#include "wsprdecoder.hh"
#include <cstdlib>
#include <cmath>
#include <map>
#include "constants.hh"
#include "fano.h"
#include <iostream>


int floatcomp(const void* elem1, const void* elem2) {
  if (*reinterpret_cast<const float*>(elem1) < *reinterpret_cast<const float*>(elem2))
    return -1;
  return *reinterpret_cast<const float*>(elem1) > *reinterpret_cast<const float*>(elem2);
}


/* ********************************************************************************************* *
 * Implementation of WSPR Decoder
 * ********************************************************************************************* */
WSPRDecoder::WSPRDecoder(MessageHandler *handler, float F0)
  : Sink(), _msgHandler(handler), _mutex(), _semaphore(), _worker(),
    _F0(F0), _Fs(1500), _N(46080)
{
  _inSync = false;
  _period = 120 * Fs;
  _input = new int16_t[_period];
  _phase = 0;

  size_t M = 32*_N;
  _inBuffer = reinterpret_cast<float *>( fftwf_malloc(sizeof(float)*M) );
  _fftOut   = reinterpret_cast<fftwf_complex*>( fftwf_malloc(sizeof(fftwf_complex)*(M/2+1)) );
  _fftIn    = reinterpret_cast<fftwf_complex*>( fftwf_malloc(sizeof(fftwf_complex)*(_N)) );
  _rfftFWD  = fftwf_plan_dft_r2c_1d(int(M), _inBuffer, _fftOut, FFTW_ESTIMATE);
  _fftREV   = fftwf_plan_dft_1d(int(_N), _fftIn, _fftOut, FFTW_BACKWARD, FFTW_ESTIMATE);
  _fftFWD   = fftwf_plan_dft_1d(512, _fftIn, _fftOut, FFTW_FORWARD, FFTW_ESTIMATE);
  _sigI   = reinterpret_cast<float *>(fftwf_malloc(sizeof(float)*_N));
  _sigQ   = reinterpret_cast<float *>(fftwf_malloc(sizeof(float)*_N));

  //static const float wdt = 2*M_PI*ModeWSPR.sym_rate/_Fs;
  for(int i=0; i<512; i++) {
    //w[i] = std::sin(wdt*i);
    w[i] = float(std::sin(0.006147931*i));
  }

  // start worker
  _workerRun = true;
  _worker = std::thread(& WSPRDecoder::_decode, this);
}

WSPRDecoder::~WSPRDecoder() {
  delete[] _input;

  fftwf_destroy_plan(_rfftFWD);
  fftwf_destroy_plan(_fftFWD);
  fftwf_destroy_plan(_fftREV);

  fftwf_free(_inBuffer);
  fftwf_free(_fftOut);
  fftwf_free(_fftIn);
  fftwf_free(_sigI);
  fftwf_free(_sigQ);

  _workerRun = false;
  _semaphore.notify_all();
  _worker.join();
}

void
WSPRDecoder::setup(float F0, bool sync, int64_t phase){
  _F0 = F0;
  _inSync = sync;
  _phase = phase;
}

void
WSPRDecoder::sync() {
  if (_inSync)
    return;

  time_t clk;
  time(&clk);
  struct tm *utc = gmtime(&clk);
  int t = (60*utc->tm_min+utc->tm_sec) % int(_period/Fs);
  _phase = int64_t(t*Fs);
  _inSync = true;
}

void
WSPRDecoder::write(const int16_t *samples, int64_t nsamples) {
  sync();

  static const int64_t len = int64_t(114*Fs);
  for (int64_t i=0; i<nsamples; i++) {
    _input[_phase++] = samples[i];
    if (len == _phase) {
      _semaphore.notify_one();
    }
    if (_period == _phase)
      _phase = 0;
  }
}

void
WSPRDecoder::_decode(WSPRDecoder *self) {
  std::cerr << "WSPR decode-worker created..." << std::endl;
  while (self->_workerRun) {
    std::unique_lock<std::mutex> lk(self->_mutex);
    self->_semaphore.wait(lk);
    if (! self->_workerRun)
      break;
    std::cerr << "WSPR decode-worker started..." << std::endl;
    self->prep_WSPR_signal();
    lk.unlock();
    self->decode();
    std::cerr << "... WSPR decode-worker finished. Wait." << std::endl;
  }
}

void
WSPRDecoder::prep_WSPR_signal() {
  size_t M = 32*_N, L=size_t(114*Fs), off=size_t(Fs*ModeWSPR.delay);
  // Skip fist second
  for (size_t i=0; i<L; i++) {
    _inBuffer[i] = float(_input[off+i])/(1<<16);
  }
  for (size_t i=L; i<M; i++) {
    _inBuffer[i] = 0;
  }
  fftwf_execute(_rfftFWD);

  double dF = Fs/M;
  size_t i0 = size_t(_F0/dF + 0.5);
  for (size_t i=0; i<_N; i++) {
    size_t j=i0+i;
    if (i>_N/2)
      j = j - _N;
    _fftIn[i][0] = _fftOut[j][0];
    _fftIn[i][1] = _fftOut[j][1];
  }
  fftwf_execute(_fftREV);

  for (size_t i=0; i<_N; i++) {
    _sigI[i] = _fftOut[i][0]/1000.0f;
    _sigQ[i] = _fftOut[i][1]/1000.0f;
  }
}

void
WSPRDecoder::decode() {
  // Do windowed ffts over 2 symbols, stepped by half symbols
  static const int nffts = 4*int(_N/512)-1, nbits = 81, symfac=50, iifac=8, delta=60;
  // First and Second sync limit
  static const float minsync1=0.10, minsync2=0.12;
  static const float df=375.0/256.0/2;
  static const float minrms=52.0 * (symfac/64.0);
  static const unsigned int maxcycles=10000;
  float fmin=-110, fmax=110;
  float dialfreq_error = 0.0;

  int nblocksize, npasses=2;
  float maxdrift;

  unsigned char symbols[2*nbits], decdata[11], channel_symbols[2*nbits];
  char callsign[13], call_loc_pow[23];
  memset(symbols, 0, 2*nbits); memset(decdata, 0, 11); memset(channel_symbols, 0, 2*nbits);
  memset(callsign, 0, 13); memset(call_loc_pow, 0, 23);
  signed char message[] = {-9,13,-35,123,57,-39,64,0,0,0,0};

  int mettab[2][256];
  // setup metric table
  for(int i=0; i<256; i++) {
      mettab[0][i]=round( 10*(metric_tables[2][i]-0.45) );
      mettab[1][i]=round( 10*(metric_tables[2][255-i]-0.45) );
  }

  std::map<std::string, int> allcalls;

  //*************** main loop starts here *****************
  for (int ipass=0; ipass<npasses; ipass++)
  {
    if (0 == ipass) {
      nblocksize = 1;
      maxdrift   = 4;
    }

    if (1 == ipass) {
      nblocksize = 3;    // try all blocksizes up to 3
      maxdrift   = 0;    // no drift for smaller frequency estimator variance
    }

    for (int i=0; i<nffts; i++) {
      for (int j=0; j<512; j++) {
        int k = i*128+j;
        if (k<int(_N)) {
          _fftIn[j][0] = _sigI[k] * w[j];
          _fftIn[j][1] = _sigQ[k] * w[j];
        } else {
          // zero pad
          _fftIn[j][0] = 0;
          _fftIn[j][1] = 0;
        }
      }
      fftwf_execute(_fftFWD);
      for (int j=0; j<512; j++) {
        int k = j+256;
        if (511 < k)
          k -= 512;
        ps[j][i] = _fftOut[k][0]*_fftOut[k][0] + _fftOut[k][1]*_fftOut[k][1];
      }
    }
    // Compute average spectrum
    for (int i=0; i<512; i++)
      psavg[i] = 0.0;
    for (int i=0; i<nffts; i++) {
      for (int j=0; j<512; j++) {
        psavg[j] += ps[j][i];
      }
    }

    // Smooth with 7-point window and limit spectrum to +/-150 Hz
    int window[7] = {1,1,1,1,1,1,1};
    float smspec[411];
    for (int i=0; i<411; i++) {
      smspec[i] = 0.0;
      for(int j=-3; j<=3; j++) {
        int k=256-205+i+j;
        smspec[i] = smspec[i] + window[j+3]*psavg[k];
      }
    }

    // Sort spectrum values, then pick off noise level as a percentile
    float tmpsort[411];
    for (int j=0; j<411; j++) {
      tmpsort[j] = smspec[j];
    }
    qsort(tmpsort, 411, sizeof(float), floatcomp);
    // Noise level of spectrum is estimated as 123/411= 30'th percentile
    float noise_level = tmpsort[122];

    // Renormalize spectrum so that (large) peaks represent an estimate of snr.
    // We know from experience that threshold snr is near -7dB in wspr bandwidth,
    // corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
    // The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15.

    // this is min snr in wspr bw
    float min_snr = std::pow(10.0f,-8.0f/10.0f), snr_scaling_factor = 26.3f;
    for (int j=0; j<411; j++) {
      smspec[j] = smspec[j]/noise_level - 1.0f;
      if( smspec[j] < min_snr)
        smspec[j] = 0.1f*min_snr;
      continue;
    }

    // Find all local maxima in smoothed spectrum.
    for (int i=0; i<200; i++) {
      freq0[i]=0.0;
      snr0[i]=0.0;
      drift0[i]=0.0;
      shift0[i]=0;
      sync0[i]=0.0;
    }

    int npk=0;
    for(int j=1; j<410; j++) {
      bool candidate = (smspec[j]>smspec[j-1]) && (smspec[j]>smspec[j+1]) && (npk<200);
      if ( candidate ) {
        freq0[npk] = (j-205)*df;
        snr0[npk] = 10*log10(smspec[j]) - snr_scaling_factor;
        npk++;
      }
    }

    // Compute corrected fmin, fmax, accounting for dial frequency error
    fmin += dialfreq_error;    // dialfreq_error is in units of Hz
    fmax += dialfreq_error;

    // Don't waste time on signals outside of the range [fmin,fmax].
    {
      int i=0;
      for(int j=0; j<npk; j++) {
        if( freq0[j] >= fmin && freq0[j] <= fmax ) {
          freq0[i]=freq0[j];
          snr0[i]=snr0[j];
          i++;
        }
      }
      npk=i;
    }

    std::cerr << " ... " << npk << " peaks in freq. range [" << fmin << "," << fmax << "]" << std::endl;
    // bubble sort on snr, bringing freq along for the ride
    for (int pass = 1; pass <= npk - 1; pass++) {
      for (int k = 0; k < npk - pass ; k++) {
        if (snr0[k] < snr0[k+1]) {
          float tmp = snr0[k];
          snr0[k] = snr0[k+1];
          snr0[k+1] = tmp;
          tmp = freq0[k];
          freq0[k] = freq0[k+1];
          freq0[k+1] = tmp;
        }
      }
    }

    // Make coarse estimates of shift (DT), freq, and drift
    // - Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
    //   to nominal start time, which is 2 seconds into the file
    //
    // - Calculates shift relative to the beginning of the file
    //
    // - Negative shifts mean that signal started before start of file
    //
    // - The program prints DT = shift-2 s
    //
    // - Shifts that cause sync vector to fall off of either end of the data
    //   vector are accommodated by "partial decoding", such that missing
    //   symbols produce a soft-decision symbol value of 128
    //
    // - The frequency drift model is linear, deviation of +/- drift/2 over the
    //   span of 162 symbols, with deviation equal to 0 at the center of the
    //   signal vector.
    int kindex;
    float p0,p1,p2,p3;
    for(int j=0; j<npk; j++) {                              //For each candidate...
      float smax = -1e30f;
      int if0 = int(freq0[j]/df+256);
      for (int ifr=if0-2; ifr<=if0+2; ifr++) {                      //Freq search
        for (int k0=-10; k0<22; k0++) {                             //Time search
          for (int idrift=-int(maxdrift); idrift<=maxdrift; idrift++) {  //Drift search
            float ss=0.0;
            float pow=0.0;
            for (int k=0; k<162; k++) {                             //Sum over symbols
              int ifd = ifr + int(float(k-81.0f)/81.0f * float(idrift)/(2.0f*df));
              kindex=k0+2*k;
              if( kindex < nffts ) {
                p0=ps[ifd-3][kindex];
                p1=ps[ifd-1][kindex];
                p2=ps[ifd+1][kindex];
                p3=ps[ifd+3][kindex];

                p0=sqrt(p0);
                p1=sqrt(p1);
                p2=sqrt(p2);
                p3=sqrt(p3);

                ss=ss+(2*sync_[k]-1)*((p1+p3)-(p0+p2));
                pow=pow+p0+p1+p2+p3;
              }
            }
            float sync1 = ss/pow;
            if (sync1 > smax) {                  //Save coarse parameters
              smax=sync1;
              shift0[j]=128*(k0+1);
              drift0[j]=idrift;
              freq0[j]=(ifr-256)*df;
              sync0[j]=sync1;
            }
          }
        }
      }
    }

    // Refine the estimates of freq, shift using sync as a metric.
    // Sync is calculated such that it is a float taking values in the range
    // [0.0,1.0].
    //
    // Function sync_and_demodulate has three modes of operation
    // mode is the last argument:
    //
    // 0 = no frequency or drift search. find best time lag.
    // 1 = no time lag or drift search. find best frequency.
    // 2 = no frequency or time lag search. Calculate soft-decision
    // symbols using passed frequency and shift.
    //
    // NB: best possibility for OpenMP may be here: several worker threads
    // could each work on one candidate at a time.
    for (int j=0; j<npk; j++) {
      memset(symbols,0,sizeof(char)*nbits*2);
      memset(callsign,0,sizeof(char)*13);
      memset(call_loc_pow,0,sizeof(char)*23);

      float f1=freq0[j], drift1=drift0[j], sync1=sync0[j];
      int shift1=shift0[j];

      // coarse-grid lag and freq search, then if sync>minsync1 continue
      float fstep=0.0, ifmin=0, ifmax=0;
      int lagmin=shift1-128;
      int lagmax=shift1+128;
      int lagstep=64;
      sync_and_demodulate(symbols, &f1, ifmin, ifmax, fstep, &shift1,
                          lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);

      fstep=0.25; ifmin=-2; ifmax=2;
      sync_and_demodulate(symbols, &f1, ifmin, ifmax, fstep, &shift1,
                          lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);

      if (ipass == 0) {
        // refine drift estimate
        fstep = 0.0; ifmin = 0; ifmax = 0;
        float driftp,driftm,syncp,syncm;
        driftp = drift1+0.5;
        sync_and_demodulate(symbols, &f1, ifmin, ifmax, fstep, &shift1,
                            lagmin, lagmax, lagstep, &driftp, symfac, &syncp, 1);

        driftm=drift1-0.5;
        sync_and_demodulate(symbols, &f1, ifmin, ifmax, fstep, &shift1,
                            lagmin, lagmax, lagstep, &driftm, symfac, &syncm, 1);

        if (syncp > sync1) {
          drift1=driftp;
          sync1=syncp;
        } else if (syncm > sync1) {
          drift1=driftm;
          sync1=syncm;
        }
      }

      // fine-grid lag and freq search
      bool worth_a_try = false;
      if( sync1 > minsync1 ) {
        lagmin = shift1-32; lagmax = shift1+32; lagstep=16;
        sync_and_demodulate(symbols, &f1, ifmin, ifmax, fstep, &shift1,
                            lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);

        // fine search over frequency
        fstep=0.05; ifmin=-2; ifmax=2;
        sync_and_demodulate(symbols, &f1, ifmin, ifmax, fstep, &shift1,
                            lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);
        worth_a_try = true;
      }

      unsigned int metric, cycles, maxnp;

      bool not_decoded=true;
      int ib=1;
      while (ib <= nblocksize && not_decoded) {
        int blocksize=ib;
        int idt=0;
        while ( worth_a_try && not_decoded && idt<=(128/iifac)) {
          int ii = (idt+1)/2;
          if (1 == idt%2)
            ii=-ii;
          ii = iifac*ii;
          int jittered_shift = shift1+ii;

          // Use mode 2 to get soft-decision symbols
          noncoherent_sequence_detection(symbols, &f1, &jittered_shift, &drift1, symfac, &blocksize);

          float sq=0.0;
          for(int i=0; i<162; i++) {
            float y = (float)symbols[i] - 128.0;
            sq += y*y;
          }
          float rms=sqrt(sq/162.0);

          if((sync1 > minsync2) && (rms > minrms)) {
            deinterleave(symbols);
            not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,nbits,
                               mettab,delta,maxcycles);
          }
          idt++;
        }
        ib++;
      }

      if (worth_a_try && !not_decoded) {
        for (int i=0; i<11; i++) {
          if (decdata[i]>127) {
            message[i]=decdata[i]-256;
          } else {
            message[i]=decdata[i];
          }
        }

        // Unpack the decoded message, update the hashtable, apply
        // sanity checks on grid and power, and return
        // call_loc_pow string and also callsign (for de-duping).
        bool noprint = unpk_(message, call_loc_pow, callsign);

        // subtract even on last pass
        if (!noprint && (ipass < npasses )) {
          if ( get_wspr_channel_symbols(call_loc_pow, channel_symbols) ) {
            subtract_signal2(f1, shift1, drift1, channel_symbols);
          } else {
            break;
          }
        }

        // Remove dupes (same callsign and freq within 3 Hz)
        bool dupe = allcalls.count(callsign) && (std::abs(f1-allcalls[callsign])<3.0);

        if (!noprint && !dupe) {
          allcalls[callsign] = f1;
          std::cerr << "Got " << call_loc_pow << std::endl;
          if (_msgHandler)
            _msgHandler->handle(call_loc_pow);
        }
      }
    }
  }
}

void
WSPRDecoder::sync_and_demodulate(unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
                                 int *shift1, int lagmin, int lagmax, int lagstep,
                                 float *drift1, int symfac, float *sync, int mode)
{
  /***********************************************************************
   * mode = 0: no frequency or drift search. find best time lag.          *
   *        1: no time lag or drift search. find best frequency.          *
   *        2: no frequency or time lag search. calculate soft-decision   *
   *           symbols using passed frequency and shift.                  *
   ************************************************************************/

  static float fplast=-10000.0;
  static float dt=1.0/375.0, df=375.0/256.0;
  static float pi=3.14159265358979323846;
  float twopidt, df15=df*1.5, df05=df*0.5;

  int i, j, k, lag;
  float i0[162],q0[162],i1[162],q1[162],i2[162],q2[162],i3[162],q3[162];
  float p0,p1,p2,p3,cmet,totp,syncmax,fac;
  float c0[256],s0[256],c1[256],s1[256],c2[256],s2[256],c3[256],s3[256];
  float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
      dphi3, cdphi3, sdphi3;
  float f0=0.0, fp, ss, fbest=0.0, fsum=0.0, f2sum=0.0, fsymb[162];
  int best_shift = 0, ifreq;

  syncmax=-1e30;
  if( mode == 0 ) {ifmin=0; ifmax=0; fstep=0.0; f0=*f1;}
  if( mode == 1 ) {lagmin=*shift1;lagmax=*shift1;f0=*f1;}
  if( mode == 2 ) {lagmin=*shift1;lagmax=*shift1;ifmin=0;ifmax=0;f0=*f1;}

  twopidt=2*pi*dt;
  for(ifreq=ifmin; ifreq<=ifmax; ifreq++) {
    f0=*f1+ifreq*fstep;
    for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) {
      ss=0.0;
      totp=0.0;
      for (i=0; i<162; i++) {
        fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
        if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
          dphi0=twopidt*(fp-df15);
          cdphi0=cos(dphi0);
          sdphi0=sin(dphi0);

          dphi1=twopidt*(fp-df05);
          cdphi1=cos(dphi1);
          sdphi1=sin(dphi1);

          dphi2=twopidt*(fp+df05);
          cdphi2=cos(dphi2);
          sdphi2=sin(dphi2);

          dphi3=twopidt*(fp+df15);
          cdphi3=cos(dphi3);
          sdphi3=sin(dphi3);

          c0[0]=1; s0[0]=0;
          c1[0]=1; s1[0]=0;
          c2[0]=1; s2[0]=0;
          c3[0]=1; s3[0]=0;

          for (j=1; j<256; j++) {
            c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
            s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
            c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
            s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
            c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
            s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
            c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
            s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
          }
          fplast = fp;
        }

        i0[i]=0.0; q0[i]=0.0;
        i1[i]=0.0; q1[i]=0.0;
        i2[i]=0.0; q2[i]=0.0;
        i3[i]=0.0; q3[i]=0.0;

        for (j=0; j<256; j++) {
          k=lag+i*256+j;
          if( (k>0) && (k<int(_N)) ) {
            i0[i]=i0[i] + _sigI[k]*c0[j] + _sigQ[k]*s0[j];
            q0[i]=q0[i] - _sigI[k]*s0[j] + _sigQ[k]*c0[j];
            i1[i]=i1[i] + _sigI[k]*c1[j] + _sigQ[k]*s1[j];
            q1[i]=q1[i] - _sigI[k]*s1[j] + _sigQ[k]*c1[j];
            i2[i]=i2[i] + _sigI[k]*c2[j] + _sigQ[k]*s2[j];
            q2[i]=q2[i] - _sigI[k]*s2[j] + _sigQ[k]*c2[j];
            i3[i]=i3[i] + _sigI[k]*c3[j] + _sigQ[k]*s3[j];
            q3[i]=q3[i] - _sigI[k]*s3[j] + _sigQ[k]*c3[j];
          }
        }
        p0 = i0[i]*i0[i] + q0[i]*q0[i];
        p1 = i1[i]*i1[i] + q1[i]*q1[i];
        p2 = i2[i]*i2[i] + q2[i]*q2[i];
        p3 = i3[i]*i3[i] + q3[i]*q3[i];

        p0=sqrt(p0);
        p1=sqrt(p1);
        p2=sqrt(p2);
        p3=sqrt(p3);

        totp=totp+p0+p1+p2+p3;
        cmet=(p1+p3)-(p0+p2);
        ss = (sync_[i] == 1) ? ss+cmet : ss-cmet;
        if( mode == 2) {                 //Compute soft symbols
          if(sync_[i]==1) {
            fsymb[i]=p3-p1;
          } else {
            fsymb[i]=p2-p0;
          }
        }
      }
      ss=ss/totp;
      if( ss > syncmax ) {          //Save best parameters
        syncmax=ss;
        best_shift=lag;
        fbest=f0;
      }
    } // lag loop
  } //freq loop

  if( mode <=1 ) {                       //Send best params back to caller
    *sync=syncmax;
    *shift1=best_shift;
    *f1=fbest;
    return;
  }

  if( mode == 2 ) {
    *sync=syncmax;
    for (i=0; i<162; i++) {              //Normalize the soft symbols
      fsum=fsum+fsymb[i]/162.0;
      f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
    }
    fac=sqrt(f2sum-fsum*fsum);
    for (i=0; i<162; i++) {
      fsymb[i]=symfac*fsymb[i]/fac;
      if( fsymb[i] > 127) fsymb[i]=127.0;
      if( fsymb[i] < -128 ) fsymb[i]=-128.0;
      symbols[i]=fsymb[i] + 128;
    }
    return;
  }
  return;
}

void
WSPRDecoder::noncoherent_sequence_detection(unsigned char *symbols, float *f1,  int *shift1,
                                            float *drift1, int symfac, int *nblocksize)
{
  /************************************************************************
   *  Noncoherent sequence detection for wspr.                            *
   *  Allowed block lengths are nblock=1,2,3,6, or 9 symbols.             *
   *  Longer block lengths require longer channel coherence time.         *
   *  The whole block is estimated at once.                               *
   *  nblock=1 corresponds to noncoherent detection of individual symbols *
   *     like the original wsprd symbol demodulator.                      *
   ************************************************************************/
  static float fplast=-10000.0;
  static float dt=1.0/375.0, df=375.0/256.0, twopidt=2*M_PI*dt;
  float df15=df*1.5, df05=df*0.5;

  int i, j, k, lag, itone, ib, b, nblock, nseq, imask;
  float xi[512],xq[512];
  float is[4][162],qs[4][162],cf[4][162],sf[4][162],cm,sm,cmp,smp;
  float p[512],fac,xm1,xm0;
  float c0[257],s0[257],c1[257],s1[257],c2[257],s2[257],c3[257],s3[257];
  float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
      dphi3, cdphi3, sdphi3;
  float f0, fp, fsum=0.0, f2sum=0.0, fsymb[162];

  f0=*f1;
  lag=*shift1;
  nblock=*nblocksize;
  nseq=1<<nblock;

  for (i=0; i<162; i++) {
    fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
    if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
      dphi0=twopidt*(fp-df15);
      cdphi0=cos(dphi0);
      sdphi0=sin(dphi0);

      dphi1=twopidt*(fp-df05);
      cdphi1=cos(dphi1);
      sdphi1=sin(dphi1);

      dphi2=twopidt*(fp+df05);
      cdphi2=cos(dphi2);
      sdphi2=sin(dphi2);

      dphi3=twopidt*(fp+df15);
      cdphi3=cos(dphi3);
      sdphi3=sin(dphi3);

      c0[0]=1; s0[0]=0;
      c1[0]=1; s1[0]=0;
      c2[0]=1; s2[0]=0;
      c3[0]=1; s3[0]=0;

      for (j=1; j<257; j++) {
        c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
        s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
        c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
        s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
        c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
        s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
        c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
        s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
      }

      fplast = fp;
    }

    cf[0][i]=c0[256]; sf[0][i]=s0[256];
    cf[1][i]=c1[256]; sf[1][i]=s1[256];
    cf[2][i]=c2[256]; sf[2][i]=s2[256];
    cf[3][i]=c3[256]; sf[3][i]=s3[256];

    is[0][i]=0.0; qs[0][i]=0.0;
    is[1][i]=0.0; qs[1][i]=0.0;
    is[2][i]=0.0; qs[2][i]=0.0;
    is[3][i]=0.0; qs[3][i]=0.0;

    for (j=0; j<256; j++) {
      k=lag+i*256+j;
      if( (k>0) && (k<int(_N)) ) {
        is[0][i] = is[0][i] + _sigI[k]*c0[j] + _sigQ[k]*s0[j];
        qs[0][i] = qs[0][i] - _sigI[k]*s0[j] + _sigQ[k]*c0[j];
        is[1][i] = is[1][i] + _sigI[k]*c1[j] + _sigQ[k]*s1[j];
        qs[1][i] = qs[1][i] - _sigI[k]*s1[j] + _sigQ[k]*c1[j];
        is[2][i] = is[2][i] + _sigI[k]*c2[j] + _sigQ[k]*s2[j];
        qs[2][i] = qs[2][i] - _sigI[k]*s2[j] + _sigQ[k]*c2[j];
        is[3][i] = is[3][i] + _sigI[k]*c3[j] + _sigQ[k]*s3[j];
        qs[3][i] = qs[3][i] - _sigI[k]*s3[j] + _sigQ[k]*c3[j];
      }
    }
  }

  for (i=0; i<162; i=i+nblock) {
    for (j=0;j<nseq;j++) {
      xi[j]=0.0; xq[j]=0.0;
      cm=1; sm=0;
      for (ib=0; ib<nblock; ib++) {
        b=(j&(1<<(nblock-1-ib)))>>(nblock-1-ib);
        itone=sync_[i+ib]+2*b;
        xi[j]=xi[j]+is[itone][i+ib]*cm + qs[itone][i+ib]*sm;
        xq[j]=xq[j]+qs[itone][i+ib]*cm - is[itone][i+ib]*sm;
        cmp=cf[itone][i+ib]*cm - sf[itone][i+ib]*sm;
        smp=sf[itone][i+ib]*cm + cf[itone][i+ib]*sm;
        cm=cmp; sm=smp;
      }
      p[j]=xi[j]*xi[j]+xq[j]*xq[j];
      p[j]=sqrt(p[j]);
    }
    for (ib=0; ib<nblock; ib++) {
      imask=1<<(nblock-1-ib);
      xm1=0.0; xm0=0.0;
      for (j=0; j<nseq; j++) {
        if((j & imask)!=0) {
          if(p[j] > xm1) xm1=p[j];
        }
        if((j & imask)==0) {
          if(p[j]>xm0) xm0=p[j];
        }
      }
      fsymb[i+ib]=xm1-xm0;
    }
  }
  for (i=0; i<162; i++) {              //Normalize the soft symbols
    fsum=fsum+fsymb[i]/162.0;
    f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
  }
  fac=std::sqrt(f2sum-fsum*fsum);
  for (i=0; i<162; i++) {
    fsymb[i]=symfac*fsymb[i]/fac;
    if( fsymb[i] > 127) fsymb[i]=127.0;
    if( fsymb[i] < -128 ) fsymb[i]=-128.0;
    symbols[i]=fsymb[i] + 128;
  }
  return;
}

void
WSPRDecoder::subtract_signal2(float f0, int shift0, float drift0, unsigned char* channel_symbols)
{
  const float dt=1.0/375.0, df=375.0/256.0, twopidt=2.0*M_PI*dt;
  const int nsym=162, nspersym=256, nfilt=256, nsig=nsym*nspersym; //nfilt must be even number.
  const int nc2=45000;

  float *refi, *refq, *ci, *cq, *cfi, *cfq;
  refi = reinterpret_cast<float *>(calloc(nc2,sizeof(float)));
  refq = reinterpret_cast<float *>(calloc(nc2,sizeof(float)));
  ci   = reinterpret_cast<float *>(calloc(nc2,sizeof(float)));
  cq   = reinterpret_cast<float *>(calloc(nc2,sizeof(float)));
  cfi  = reinterpret_cast<float *>(calloc(nc2,sizeof(float)));
  cfq  = reinterpret_cast<float *>(calloc(nc2,sizeof(float)));

  /* ***************************************************************************** *
   * Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
   * Reference is:                       r(t) = exp( j*phi(t) )
   * Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
   * so c(t) has phase angle theta-phi
   * Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
   * ***************************************************************************** */

    // create reference wspr signal vector, centered on f0.
  for (int i=0; i<nsym; i++) {
    float cs=(float)channel_symbols[i];

    float dphi = twopidt * (
          f0 + (drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0)
          + (cs-1.5)*df );

    float phi=0;
    for (int j=0; j<nspersym; j++ ) {
      int ii = nspersym*i+j;
      refi[ii] = std::cos(phi); //cannot precompute sin/cos because dphi is changing
      refq[ii] = std::sin(phi);
      phi += dphi;
    }
  }

  // s(t) * conjugate(r(t))
  // beginning of first symbol in reference signal is at i=0
  // beginning of first symbol in received data is at shift0.
  // filter transient lasts nfilt samples
  // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
  for (int i=0; i<nsym*nspersym; i++) {
    int k = shift0+i;
    if( (k>0) && (k<int(_N)) ) {
      ci[i+nfilt] = _sigI[k]*refi[i] + _sigQ[k]*refq[i];
      cq[i+nfilt] = _sigQ[k]*refi[i] - _sigI[k]*refq[i];
    }
  }

  //lowpass filter and remove startup transient
  float w[nfilt], norm=0, partialsum[nfilt];
  for (int i=0; i<nfilt; i++)
    partialsum[i]=0.0;
  for (int i=0; i<nfilt; i++) {
    w[i] = std::sin(M_PI*(float)i/(float)(nfilt-1));
    norm=norm+w[i];
  }
  for (int i=0; i<nfilt; i++) {
    w[i]=w[i]/norm;
  }
  for (int i=1; i<nfilt; i++) {
    partialsum[i]=partialsum[i-1]+w[i];
  }

  // LPF
  for (int i=nfilt/2; i<45000-nfilt/2; i++) {
    cfi[i]=0.0; cfq[i]=0.0;
    for (int j=0; j<nfilt; j++) {
      cfi[i]=cfi[i]+w[j]*ci[i-nfilt/2+j];
      cfq[i]=cfq[i]+w[j]*cq[i-nfilt/2+j];
    }
  }

  // subtract c(t)*r(t) here
  // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq)+cq*refi)
  // beginning of first symbol in reference signal is at i=nfilt
  // beginning of first symbol in received data is at shift0.
  for (int i=0; i<nsig; i++) {
    if( i<nfilt/2 ) {        // take care of the end effect (LPF step response) here
      norm=partialsum[nfilt/2+i];
    } else if( i>(nsig-1-nfilt/2) ) {
      norm=partialsum[nfilt/2+nsig-1-i];
    } else {
      norm=1.0;
    }
    int k=shift0+i;
    int j=i+nfilt;
    if( (k>0) && (k<int(_N)) ) {
      _sigI[k] -= (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
      _sigQ[k] -= (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
    }
  }

  free(refi);
  free(refq);
  free(ci);
  free(cq);
  free(cfi);
  free(cfq);

  return;
}
