#ifndef CWSPR_TRAITS_H
#define CWSPR_TRAITS_H

#define Fs 12000.


typedef enum {
	MODE_JT4, MODE_WSPR
} ModeId;

typedef struct {
	ModeId mode;
	int symbols;
	double period;
	double delay;
	double sym_rate;
	double fsk_shift;
} Mode;


#endif // CWSPR_TRAITS_H
