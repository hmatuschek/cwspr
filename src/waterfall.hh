#ifndef WATERFALL_HH
#define WATERFALL_HH

#include <QLabel>
#include "interfaces.hh"
#include <fftw3.h>
#include <QTimer>
#include <QImage>

class Application;

class WaterfallSink: public QObject, public Sink
{
	Q_OBJECT

public:
	explicit WaterfallSink(int linelen, QObject *parent=nullptr);
	virtual ~WaterfallSink();

	void write(const int16_t *samples, qint64 nsamples);

	int nFFT() const;
	const double *spec() const;

protected:
	int _Ns, _lineLen, _Nf, _nAvg;
	fftwf_complex *_sigBuff;
	double *_spec, *_aspec;
	fftwf_plan _plan;
};


class Waterfall : public QLabel
{
  Q_OBJECT

public:
	explicit Waterfall(Application &app, QWidget *parent = nullptr);

	WaterfallSink *sink();

protected slots:
	void update();

protected:
	void resizeEvent(QResizeEvent *event);

protected:
	WaterfallSink *_sink;
	QImage _plot;
	QTimer _updateTimer;
};

#endif // WATERFALL_HH
