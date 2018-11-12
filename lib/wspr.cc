#include "wspr.hh"
#include <iostream>
#include "traits.hh"
#include "fano.h"
#include <cmath>


extern "C" {
void encode232_(const char *data, int *dlen, char *sym);
void inter_wspr_(char *sym, int *dir);
void wqencode_(const char *msg, int *ntype, char *data);
}


const int
WSPR::sync_[162] = {
  1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
  0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
  0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
  1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
  0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
  0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
  0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
  0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
  0,0 };

Mode ModeWSPR = { MODE_WSPR, 162, 120.0, 1.0, 1.464738, 1.4648};

int
WSPR::gen_wspr(const std::string &message, std::vector<int> &symbols) {
  char msg[22], data[13];
  memset(msg, ' ', 22);
  memcpy(msg, message.c_str(), std::min(size_t(22), message.size()));
  memset(data, 0, 13);

  int ntype=0;
  wqencode_(msg, &ntype, data);

  char sym[206];
  memset(sym, 0, 206);
  int slen=206, dir=1;
  encode232_(data, &slen, sym);
  inter_wspr_(sym, &dir);

  symbols.resize(162);
  for (int i=0; i<162; i++) {
    symbols[i] = 2*int(sym[i])+sync_[i];
  }

  return ntype;
}


void
WSPR::interleave(unsigned char *sym)
{
  unsigned char tmp[162];
  unsigned char p=0, i=0, j;

  while (p<162) {
    j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
    if (j < 162 ) {
      tmp[j]=sym[p];
      p++;
    }
    i++;;
  }
  for (i=0; i<162; i++) {
    sym[i]=tmp[i];
  }
}

void
WSPR::deinterleave(unsigned char *sym)
{
  unsigned char tmp[162];
  unsigned char p=0, i=0, j;

  while (p<162) {
    j = ((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
    if (j < 162 ) {
      tmp[p]=sym[j];
      p++;
    }
    i++;
  }
  for (i=0; i<162; i++) {
    sym[i]=tmp[i];
  }
}

bool
WSPR::unpk_(signed char *message, char *call_loc_pow, char *callsign)
{
  int n1,n2,ndbm;
  bool noprint=false;
  char grid[5], cdbm[3];

  unpack50(message,&n1,&n2);
  if (!unpackcall(n1,callsign))
    return true;
  if (!unpackgrid(n2, grid))
    return true;
  int ntype = (n2&127) - 64;
  callsign[12]=0;
  grid[4]=0;

  /* Based on the value of ntype, decide whether this is a Type 1, 2, or
   * 3 message.
   *
   * - Type 1: 6 digit call, grid, power - ntype is positive and is a member
   *   of the set {0,3,7,10,13,17,20...60}
   *
   * - Type 2: extended callsign, power - ntype is positive but not
   *   a member of the set of allowed powers
   *
   * - Type 3: hash, 6 digit grid, power - ntype is negative. */
  if( (ntype >= 0) && (ntype <= 62) ) {
    int nu=ntype%10;
    if( nu == 0 || nu == 3 || nu == 7 ) {
      ndbm=ntype;
      memset(call_loc_pow,0,sizeof(char)*23);
      sprintf(cdbm,"%2d",ndbm);
      strncat(call_loc_pow,callsign,strlen(callsign));
      strncat(call_loc_pow," ",1);
      strncat(call_loc_pow,grid,4);
      strncat(call_loc_pow," ",1);
      strncat(call_loc_pow,cdbm,2);
      strncat(call_loc_pow,"\0",1);
    } else {
      return true;
    }
  } else if ( ntype == -64 ) {
    // I don't know what to do with these... They show up as "A000AA" grids.
    return true;
  }

  return noprint;
}

int
WSPR::get_wspr_channel_symbols(char* rawmessage, unsigned char* symbols)
{
  int m=0, ntype=0;
  long unsigned int n=0;
  int i;
  int nu[10]={0,-1,1,0,-1,2,1,0,-1,1};
  char *callsign, *grid, *powstr;
  char grid4[5], message[23];

  memset(message,0,sizeof(char)*23);
  i=0;
  while ( rawmessage[i] != 0 && i<23 ) {
    message[i]=rawmessage[i];
    i++;
  }

  size_t i1=strcspn(message," ");
  size_t i2=strcspn(message,"/");
  size_t i3=strcspn(message,"<");
  size_t mlen=strlen(message);

  // Use the presence and/or absence of "<" and "/" to decide what
  // type of message. No sanity checks! Beware!

  if( i1 > 3 && i1 < 7 && i2 == mlen && i3 == mlen ) {
    // Type 1 message: K9AN EN50 33
    //                 xxnxxxx xxnn nn
    callsign = strtok(message," ");
    grid = strtok(NULL," ");
    powstr = strtok(NULL," ");
    int power = atoi(powstr);
    n = pack_call(callsign);
    for (i=0; i<4; i++) {
      grid4[i] = get_locator_character_code(*(grid+i));
    }
    m = pack_grid4_power(grid4,power);

  } else if ( i2 < mlen ) {  // just looks for a right slash
    // Type 2: PJ4/K1ABC 37
    callsign = strtok (message," ");
    if( i2==0 || i2>strlen(callsign) ) return 0; //guards against pathological case
    powstr = strtok (NULL," ");
    int power = atoi (powstr);
    if( power < 0 ) power=0;
    if( power > 60 ) power=60;
    power=power+nu[power%10];
    int n1, ng, nadd;
    pack_prefix(callsign, &n1, &ng, &nadd);
    ntype=power + 1 + nadd;
    m=128*ng+ntype+64;
    n=n1;
  } else {
    return 0;
  }

  // pack 50 bits + 31 (0) tail bits into 11 bytes
  unsigned char it, data[11];
  memset(data,0,sizeof(char)*11);
  it=0xFF & (n>>20);
  data[0]=it;
  it=0xFF & (n>>12);
  data[1]=it;
  it=0xFF & (n>>4);
  data[2]=it;
  it= ((n&(0x0F))<<4) + ((m>>18)&(0x0F));
  data[3]=it;
  it=0xFF & (m>>10);
  data[4]=it;
  it=0xFF & (m>>2);
  data[5]=it;
  it=(m & 0x03)<<6 ;
  data[6]=it;
  data[7]=0;
  data[8]=0;
  data[9]=0;
  data[10]=0;

  // make sure that the 11-byte data vector is unpackable
  // unpack it with the routine that the decoder will use and display
  // the result. let the operator decide whether it worked.

  char check_call_loc_pow[23], check_callsign[13];
  signed char check_data[11];
  memcpy(check_data,data,sizeof(char)*11);

  unpk_(check_data,check_call_loc_pow,check_callsign);

  unsigned int nbytes=11; // The message with tail is packed into almost 11 bytes.
  unsigned char channelbits[nbytes*8*2]; /* 162 rounded up */
  memset(channelbits,0,sizeof(char)*nbytes*8*2);

  encode(channelbits,data,nbytes);

  interleave(channelbits);

  for (i=0; i<162; i++) {
    symbols[i] = 2*channelbits[i]+sync_[i];
  }
  return 1;
}



long unsigned int
WSPR::pack_grid4_power(char const *grid4, int power) {
  long unsigned int m = (179-10*grid4[0]-grid4[2])*180 + 10*grid4[1] + grid4[3];
  m=m*128+power+64;
  return m;
}

void
WSPR::unpack50(signed char *dat, int32_t *n1, int32_t *n2)
{
  int32_t i,i4;

  i=dat[0];
  i4=i&255;
  *n1=i4<<20;

  i=dat[1];
  i4=i&255;
  *n1=*n1+(i4<<12);

  i=dat[2];
  i4=i&255;
  *n1=*n1+(i4<<4);

  i=dat[3];
  i4=i&255;
  *n1=*n1+((i4>>4)&15);
  *n2=(i4&15)<<18;

  i=dat[4];
  i4=i&255;
  *n2=*n2+(i4<<10);

  i=dat[5];
  i4=i&255;
  *n2=*n2+(i4<<2);

  i=dat[6];
  i4=i&255;
  *n2=*n2+((i4>>6)&3);
}

int
WSPR::unpackcall( int32_t ncall, char *call )
{
  char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
            'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
            'U','V','W','X','Y','Z',' '};
  int32_t n;
  int i;
  char tmp[7];

  n=ncall;
  strcpy(call,"......");
  if (n < 262177560 ) {
    i=n%27+10;
    tmp[5]=c[i];
    n=n/27;
    i=n%27+10;
    tmp[4]=c[i];
    n=n/27;
    i=n%27+10;
    tmp[3]=c[i];
    n=n/27;
    i=n%10;
    tmp[2]=c[i];
    n=n/10;
    i=n%36;
    tmp[1]=c[i];
    n=n/36;
    i=n;
    tmp[0]=c[i];
    tmp[6]='\0';
    // remove leading whitespace
    for(i=0; i<5; i++) {
      if( tmp[i] != c[36] )
        break;
    }
    sprintf(call,"%-6s",&tmp[i]);
    // remove trailing whitespace
    for(i=0; i<6; i++) {
      if( call[i] == c[36] ) {
        call[i]='\0';
      }
    }
  } else {
    return 0;
  }
  return 1;
}

int
WSPR::unpackgrid(int32_t ngrid, char *grid)
{
  char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
            'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
            'U','V','W','X','Y','Z',' '};
  int dlat, dlong;

  ngrid=ngrid>>7;
  if( ngrid < 32400 ) {
    dlat=(ngrid%180)-90;
    dlong=(ngrid/180)*2 - 180 + 2;
    if( dlong < -180 )
      dlong=dlong+360;
    if( dlong > 180 )
      dlong=dlong+360;
    int nlong = 60.0*(180.0-dlong)/5.0;
    int n1 = nlong/240;
    int n2 = (nlong - 240*n1)/24;
    grid[0] = c[10+n1];
    grid[2]=  c[n2];

    int nlat = 60.0*(dlat+90)/2.5;
    n1 = nlat/240;
    n2 = (nlat-240*n1)/24;
    grid[1]=c[10+n1];
    grid[3]=c[n2];
  } else {
    strcpy(grid,"XXXX");
    return 0;
  }
  return 1;
}

int
WSPR::unpackpfx(int32_t nprefix, char *call)
{
  char nc, pfx[4]={'\0'}, tmpcall[7];
  int i;
  int32_t n;

  strcpy(tmpcall,call);
  if( nprefix < 60000 ) {
    // add a prefix of 1 to 3 characters
    n=nprefix;
    for (i=2; i>=0; i--) {
      nc=n%37;
      if( (nc >= 0) & (nc <= 9) ) {
        pfx[i]=nc+48;
      }
      else if( (nc >= 10) & (nc <= 35) ) {
        pfx[i]=nc+55;
      }
      else {
        pfx[i]=' ';
      }
      n=n/37;
    }

    char * p = strrchr(pfx,' ');
    strcpy(call, p ? p + 1 : pfx);
    strncat(call,"/",1);
    strncat(call,tmpcall,strlen(tmpcall));

  } else {
    // add a suffix of 1 or 2 characters
    nc=nprefix-60000;
    if( (nc >= 0) & (nc <= 9) ) {
      pfx[0]=nc+48;
      strcpy(call,tmpcall);
      strncat(call,"/",1);
      strncat(call,pfx,1);
    }
    else if( (nc >= 10) & (nc <= 35) ) {
      pfx[0]=nc+55;
      strcpy(call,tmpcall);
      strncat(call,"/",1);
      strncat(call,pfx,1);
    }
    else if( (nc >= 36) & (nc <= 125) ) {
      pfx[0]=(nc-26)/10+48;
      pfx[1]=(nc-26)%10+48;
      strcpy(call,tmpcall);
      strncat(call,"/",1);
      strncat(call,pfx,2);
    }
    else {
      return 0;
    }
  }
  return 1;
}

long unsigned int
WSPR::pack_call(char const *callsign) {
  unsigned int i;
  long unsigned int n;
  char call6[6];
  memset(call6,' ',sizeof(call6));
  // callsign is 6 characters in length. Exactly.
  size_t call_len = strlen(callsign);
  if( call_len > 6 ) {
    return 0;
  }
  if( isdigit(callsign[2]) ) {
    for (i=0; i<call_len; i++) {
      call6[i]=callsign[i];
    }
  } else if( isdigit(callsign[1]) ) {
    for (i=1; i<call_len+1; i++) {
      call6[i]=callsign[i-1];
    }
  }
  for (i=0; i<6; i++) {
    call6[i]=get_callsign_character_code(call6[i]);
  }
  n = call6[0];
  n = n*36+call6[1];
  n = n*10+call6[2];
  n = n*27+call6[3]-10;
  n = n*27+call6[4]-10;
  n = n*27+call6[5]-10;
  return n;
}

char
WSPR::get_locator_character_code(char ch) {
  if( ch >=48 && ch <=57 ) { //0-9
    return ch-48;
  }
  if( ch == 32 ) {  //space
    return 36;
  }
  if( ch >= 65 && ch <= 82 ) { //A-Z
    return ch-65;
  }
  return -1;
}

char
WSPR::get_callsign_character_code(char ch) {
  if( ch >=48 && ch <=57 ) { //0-9
    return ch-48;
  }
  if( ch == 32 ) {  //space
    return 36;
  }
  if( ch >= 65 && ch <= 90 ) { //A-Z
    return ch-55;
  }
  return -1;
}

void
WSPR::pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd ) {
  size_t i;
  char call6[7]; memset(call6, 0, 7);
  size_t i1 = strcspn(callsign,"/");

  if( callsign[i1+2] == 0 ) {
    //single char suffix
    for (i=0; i<i1; i++) {
      call6[i]=callsign[i];
    }
    call6[i] = '\0';
    *n=pack_call(call6);
    *nadd=1;
    int nc = callsign[i1+1];
    if( nc >= 48 && nc <= 57 ) {
      *m=nc-48;
    } else if ( nc >= 65 && nc <= 90 ) {
      *m=nc-65+10;
    } else {
      *m=38;
    }
    *m=60000-32768+*m;
  } else if( callsign[i1+3]==0 ) {
    //two char suffix
    for (i=0; i<i1; i++) {
      call6[i]=callsign[i];
    }
    *n=pack_call(call6);
    *nadd=1;
    *m=10*(callsign[i1+1]-48)+(callsign[i1+2]-48);
    *m=60000 + 26 + *m;
  } else {
    char const * pfx = strtok (callsign,"/");
    char const * call = strtok(NULL," ");
    *n = pack_call (call);
    size_t plen=strlen (pfx);
    if( plen ==1 ) {
      *m=36;
      *m=37*(*m)+36;
    } else if( plen == 2 ) {
      *m=36;
    } else {
      *m=0;
    }
    for (i=0; i<plen; i++) {
      int nc = callsign[i];
      if( nc >= 48 && nc <= 57 ) {
        nc=nc-48;
      } else if ( nc >= 65 && nc <= 90 ) {
        nc=nc-65+10;
      } else {
        nc=36;
      }
      *m=37*(*m)+nc;
    }
    *nadd=0;
    if( *m > 32768 ) {
      *m=*m-32768;
      *nadd=1;
    }
  }
}


