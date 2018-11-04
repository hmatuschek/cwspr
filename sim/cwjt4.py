from struct import unpack,calcsize
from numpy import *
from scipy.signal import resample_poly
import scipy.io.wavfile

# Free text
stext = ("2 2 2 1 1 0 0 2 1 3 2 1 1 0 2 3 2 3 0 0 0 0 0 0 0 3 1 2 0 0 0 0 2 0 0 2 2 0 0 3 0 1 1 0 1 "
         "1 0 1 2 3 3 3 3 1 0 3 0 2 0 1 0 2 3 0 0 1 3 1 3 3 0 2 0 1 0 3 2 0 0 1 1 3 1 2 1 1 2 2 3 2 "
		 "2 0 1 3 2 1 0 1 0 1 2 3 2 1 1 3 3 1 0 1 2 1 2 1 3 2 1 2 3 0 3 3 3 2 0 1 2 1 3 0 1 3 3 1 2 "
		 "0 0 0 3 1 0 3 1 2 0 2 1 3 3 0 3 3 3 0 3 3 3 0 2 1 0 0 0 1 1 2 1 3 0 2 1 2 0 0 3 1 1 1 3 1 "
		 "0 0 3 1 2 0 0 0 1 1 0 0 2 1 2 1 1 2 3 3 3 3 2 3 0 1 0");

symbols = list(map(int, stext.split(" ")))
print(len(symbols))
ftab = 600.+array([0.,1.,2.,3.])*4.375


Fs = 12000.;
T  = 46.
dur = 1./4.375
Ns = int(Fs*T);
t  = linspace(0,T, Ns)
dt = T/Ns
Noff = int(1*Fs)
wpm = 18;
ditlen = 60./50/wpm

code = ("-.. .-.. --... -..- -.-- --..   -.. .   -.. -- ...-- -- .- -   "
		"-...-   - -. -..-   ..-. . .-.   -.-. .- .-.. .-..   -.. .-.   --- --   "
		"-...-   ..- .-.   .-. ... -  ..... -. -.   ..... -. -.   "
		"-...-   -. .- -- .   .... .- -. -. . ...   .... .- -. -. . ...   "
		"-...-   --.- - ....   -. .-.   -... . .-. .-.. .. -.   -. .-.   -... . .-. .-.. .. -.   "
		"-...-  .... .-- ..--..   "
		".-.-.   -.. -- ...-- -- .- -   -.. .   -.. -- ...-- -- .- -   -...-.-   ");

cw = code.replace("-","1110").replace(".","10").replace(" ","00")

#sig = zeros((Ns+int(9.4*Fs),))
sig = random.randn(Ns+int((60-T)*Fs))

for i in range(Ns):
	idx = symbols[int(t[i]/dur)]
	cwi = int(t[i]/ditlen)
	if ((cwi//len(cw)) % 2):
		cwi = len(cw)-1;
	else:
		cwi = cwi % len(cw)
	if ('1' == cw[cwi]):
		sig[Noff+i] += sin(2*pi*ftab[idx]*t[i])/2

sig = array(7e3*sig, dtype="int16")


scipy.io.wavfile.write("cwjt4_qso-3db.wav", Fs, sig);
