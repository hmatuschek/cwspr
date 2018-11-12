#ifndef CWSPR_WSPR_HH
#define CWSPR_WSPR_HH

#include <string>
#include <vector>
#include "traits.hh"
#include <fftw3.h>
#include <list>
#include <string>


class WSPR
{
public:
	static int get_wspr_channel_symbols(char* rawmessage, unsigned char* symbols);
	static void interleave(unsigned char *sym);
	static void pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd );
	static long unsigned int pack_call(char const *callsign);
	static long unsigned int pack_grid4_power(char const *grid4, int power);

	static void deinterleave(unsigned char *sym);
	static void unpack50(signed char *dat, int32_t *n1, int32_t *n2);
	static int unpackpfx(int32_t nprefix, char *call);
	static int unpackcall(int32_t ncall, char *call);
	static int unpackgrid(int32_t ngrid, char *grid);
	static bool unpk_(signed char *message, char *call_loc_pow, char *callsign);

	static char get_locator_character_code(char ch);
	static char get_callsign_character_code(char ch);

	static int gen_wspr(const std::string &message, std::vector<int> &symbols);

protected:
	static const int sync_[162];
};


extern Mode ModeWSPR;

#endif // CWSPR_WSPR_HH
