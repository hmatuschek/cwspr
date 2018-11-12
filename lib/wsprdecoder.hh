#ifndef WSPRDECODER_HH
#define WSPRDECODER_HH

#include "interfaces.hh"
#include "wspr.hh"
#include <condition_variable>
#include <semaphore.h>
#include <thread>


class WSPRDecoder: public WSPR, public Sink
{
public:
	WSPRDecoder(MessageHandler *handler, float F0);
	virtual ~WSPRDecoder();

	void setup(float F0, bool sync=false, int64_t phase=0);
	void write(const int16_t *samples, int64_t nsamples);

protected:
	void sync();

	void prep_WSPR_signal();
	void decode();
	void sync_and_demodulate(unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
	                         int *shift1, int lagmin, int lagmax, int lagstep,
	                         float *drift1, int symfac, float *sync, int mode);
	void noncoherent_sequence_detection(unsigned char *symbols, float *f1,  int *shift1,
	                                    float *drift1, int symfac, int *nblocksize);
	void subtract_signal2(float f0, int shift0, float drift0, unsigned char* channel_symbols);

	static void _decode(WSPRDecoder *self);

protected:
  MessageHandler *_msgHandler;
	bool _inSync;
	int16_t *_input;
	int64_t _phase;
	int64_t _period;

	std::mutex _mutex;
	std::condition_variable _semaphore;
  std::thread _worker;
  bool _workerRun;

	float _F0, _Fs;
	const size_t _N;
	fftwf_plan _rfftFWD, _fftFWD, _fftREV;
	float *_inBuffer;
	fftwf_complex *_fftIn, *_fftOut;
	float *_sigI, *_sigQ;
  float ps[512][359];
	std::list<std::string> _messages;
};

#endif // WSPRDECODER_HH
