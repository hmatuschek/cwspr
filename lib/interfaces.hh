#ifndef INTERFACES_HH
#define INTERFACES_HH

#include <cinttypes>
#include <string>


class MessageHandler
{
protected:
	explicit MessageHandler();

public:
	virtual ~MessageHandler();

	virtual void handle(const std::string &message) = 0;
};


class Sink
{
protected:
	explicit Sink();

public:
	virtual ~Sink();

	virtual void write(const int16_t *samples, int64_t nsample) = 0;
};


class Source
{
protected:
	explicit Source();

public:
	virtual ~Source();

	virtual void read(int16_t *samples, int64_t nsamples) = 0;
};


#endif // INTERFACES_HH
