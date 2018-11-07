#ifndef SOURCE_HH
#define SOURCE_HH

#include <QObject>

class Source : public QObject
{
  Q_OBJECT

public:
	explicit Source(QObject *parent = nullptr);
	virtual ~Source();

	virtual void read(int16_t *samples, qint64 nsamples) = 0;
};

#endif // SOURCE_HH
