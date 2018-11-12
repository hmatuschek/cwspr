#include "waterfall.hh"
#include "application.hh"
#include <QDebug>
#include <cmath>


/* ********************************************************************************************* *
 * Implementation of WaterfallSink
 * ********************************************************************************************* */
WaterfallSink::WaterfallSink(int lineLen, QObject *parent)
  : QObject(parent), Sink(), _Ns(0), _lineLen(lineLen), _Nf(0), _nAvg(Fs/_lineLen/2)
{
  _sigBuff = (fftwf_complex *)(fftwf_malloc(2*_lineLen*sizeof(fftwf_complex)));
  _aspec   = new double[_lineLen];
  _spec    = new double[_lineLen];
  _plan = fftwf_plan_dft_1d(2*_lineLen, _sigBuff, _sigBuff, FFTW_FORWARD, FFTW_ESTIMATE);
}

WaterfallSink::~WaterfallSink() {
  fftwf_free(_sigBuff);
  delete [] _spec;
  delete [] _aspec;
  fftwf_destroy_plan(_plan);
}

int
WaterfallSink::nFFT() const {
  return _lineLen;
}

const double *
WaterfallSink::spec() const {
  return _spec;
}

void
WaterfallSink::write(const int16_t *samples, qint64 nsamples)
{
  for (qint64 i=0; i<nsamples; i++) {
    _sigBuff[_Ns][0] = float(samples[i])/(1<<15);
    _sigBuff[_Ns][1] = 0;
    _Ns++;

    if ((2*_lineLen) == _Ns) {
      _Ns = 0;
      fftwf_execute(_plan);

      for (int j=0; j<_lineLen; j++) {
        _aspec[j] += (_sigBuff[j][0]*_sigBuff[j][0] + _sigBuff[j][1]*_sigBuff[j][1])/_lineLen/_lineLen;
      }
      _Nf++;

      if (_Nf == _nAvg) {
        _Nf = 0;
        for (int j=0; j<_lineLen; j++) {
          _spec[j] = _aspec[j]/_nAvg;
          _aspec[j] = 0;
        }
      }
    }
  }
}



/* ********************************************************************************************* *
 * Implementation of WaterfallSink
 * ********************************************************************************************* */
uint32_t colormap(double db) {
  uint32_t cols[6] = { 0xff000000, 0xff0000ff, 0xffff00ff, 0xffffffff};
  double min = -80, max=0;
  db = std::min(max, std::max(min, db));
  db -= min; db /= (max-min)/3;
  int idx = int(db);
  double frac = db-idx;
  return cols[idx] + uint32_t((cols[idx+1]-cols[idx])*frac);
}


Waterfall::Waterfall(Application &app, QWidget *parent)
    : QLabel(parent), _sink(app.waterfall()),
      _plot(512, 150, QImage::Format_RGB32), _updateTimer()
{
  setMinimumSize(512, 150);
  _plot.fill(0);
  setScaledContents(true);

  _updateTimer.setInterval(1000);
  _updateTimer.setSingleShot(false);
  connect(&_updateTimer, SIGNAL(timeout()), this, SLOT(update()));
  _updateTimer.start();

  update();
}

WaterfallSink *
Waterfall::sink() {
  return _sink;
}

void
Waterfall::update() {
  double spec[_sink->nFFT()];
  memcpy(spec, _sink->spec(), _sink->nFFT()*sizeof(double));
  memcpy(_plot.scanLine(0), _plot.scanLine(1), (_plot.height()-1)*_plot.bytesPerLine());

  double F0 = 300, dF = (900-300)/_plot.width(), F;
  for (int i=0; i<_plot.width(); i++) {
    F = F0 + i*dF;
    int idx = int(2*_sink->nFFT()*F/Fs);
    double db = 10*std::log10(spec[idx]);
    _plot.setPixel(i,(_plot.height()-1),colormap(db));
  }
  this->setPixmap(QPixmap::fromImage(_plot));
}

void
Waterfall::resizeEvent(QResizeEvent *evt) {
  _plot = QImage(512, 150, QImage::Format_RGB32);
  _plot.fill(0);
}
