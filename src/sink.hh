#ifndef SINK_HH
#define SINK_HH

#include <QObject>

class Sink: public QObject
{
	Q_OBJECT

public:
	explicit Sink(QObject *parent=nullptr);
	virtual ~Sink();

	virtual void write(const int16_t *samples, qint64 nsample) = 0;
};

#endif // SINK_HH
