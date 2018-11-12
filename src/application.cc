#include "application.hh"
#include "jt4.hh"
#include "wspr.hh"
#include <QDebug>
#include <QSettings>

Application::Application(int &argc, char *argv[])
  : QApplication(argc, argv), _encoder(nullptr), _keyer(nullptr), _waterfall(nullptr), _ptt(false),
    _pttTimer()
{
  Pa_Initialize();

  setOrganizationDomain("io.github.hmatuschek");
  setApplicationName("cwspr");

  qDebug() << "Use PortAudio" << Pa_GetVersionInfo();
  PaStreamParameters fmt = {
    Pa_GetDefaultOutputDevice(), // Output device
    1,                           // channels,
    paInt16,                     // type,
    16e-3,                       // latency
    nullptr
  };

  if (Pa_IsFormatSupported(nullptr, &fmt, Fs))
    qDebug() << "Warning: Audio format is not supported by device"
             << Pa_GetDeviceInfo(fmt.device)->name;

  if (int err = Pa_OpenStream(&_out, nullptr, &fmt, Fs, 2, paPrimeOutputBuffersUsingStreamCallback,
                              &Application::_pa_outstream_callback, this)) {
    qDebug() << "Error: Audio stream cannot be created:"
             << Pa_GetErrorText(err);
    return;
  }


  fmt = {
    Pa_GetDefaultInputDevice(),  // Input device
    1,                           // channels,
    paInt16,                     // type,
    64e-3,                       // latency
    nullptr
  };

  if (Pa_IsFormatSupported(nullptr, &fmt, Fs))
    qDebug() << "Warning: Audio format is not supported by device"
             << Pa_GetDeviceInfo(fmt.device)->name;

  if (int err = Pa_OpenStream(&_in, &fmt, nullptr, Fs, 1, paPrimeOutputBuffersUsingStreamCallback,
                              &Application::_pa_instream_callback, this)) {
    qDebug() << "Error: Audio stream cannot be created:"
             << Pa_GetErrorText(err);
    return;
  }

  QSettings settings;
  ModeId mode = ModeId(settings.value("mode", MODE_WSPR).toUInt());
  QString message = settings.value("message", "DM3MAT JO62 37").toString();
  double freq = settings.value("freq", 600.0).toDouble();

  _encoder = new Encoder(mode, message, freq, this);
  _decoder = new WSPRDecoder(nullptr, freq);
  _keyer = new Keyer(_encoder, this);
  _waterfall = new WaterfallSink(2048, this);
  Pa_StartStream(_in);

  _pttTimer.setInterval(500);
  _pttTimer.setSingleShot(true);
  connect(&_pttTimer, SIGNAL(timeout()), this, SLOT(onPTTTimeOut()));
}

Application::~Application() {
  Pa_StopStream(_in);
  Pa_Terminate();
}

void
Application::info() {
}

WaterfallSink *
Application::waterfall() {
  return _waterfall;
}

void Application::start() {
  if (Pa_IsStreamActive(_out))
    return;

  QSettings settings;
  ModeId mode = ModeId(settings.value("mode", MODE_WSPR).toUInt());
  QString call = settings.value("call", "DM3MAT").toString();
  QString loc = settings.value("locator", "JO62").toString();
  int dBm = settings.value("power", 37).toInt();
  QString message = settings.value("message", "DM3MAT JO62 37").toString();
  if (MODE_WSPR == mode)
    message = QString("%1 %2 %3").arg(call).arg(loc).arg(dBm);
  double freq = settings.value("freq", 600.0).toDouble();

  _encoder->setup(mode, message, freq);
  Pa_StartStream(_out);
}

void Application::stop() {
  Pa_StopStream(_out);
}

void
Application::key(bool down) {
  _keyer->key(down);

  if (down) {
    if (! _ptt) {
      emit ptt(true);
    }
    _ptt = true;
    _pttTimer.start();
  }
}

void
Application::ditKey(bool down) {
  _keyer->ditKey(down);

  if (down) {
    if (! _ptt) {
      emit ptt(true);
    }
    _ptt = true;
    _pttTimer.start();
  }
}

void
Application::daKey(bool down) {
  _keyer->daKey(down);

  if (down) {
    if (! _ptt) {
      emit ptt(true);
    }
    _ptt = true;
    _pttTimer.start();
  }
}

void
Application::onPTTTimeOut() {
  _ptt = false;
  emit ptt(false);
}

int
Application::_pa_outstream_callback(const void *in, void *out, unsigned long frameCount,
                                    const PaStreamCallbackTimeInfo *tInfo,
                                    PaStreamCallbackFlags status, void *userData)
{
  Q_UNUSED(in); Q_UNUSED(tInfo); Q_UNUSED(status);
  reinterpret_cast<Application *>(userData)
      ->_keyer->read(static_cast<int16_t *>(out), qint64(frameCount));
  return paContinue;
}

int
Application::_pa_instream_callback(const void *in, void *out, unsigned long frameCount,
                                   const PaStreamCallbackTimeInfo *tInfo,
                                   PaStreamCallbackFlags status, void *userData)
{
  Q_UNUSED(out); Q_UNUSED(tInfo); Q_UNUSED(status);
  reinterpret_cast<Application *>(userData)
      ->_waterfall->write(static_cast<const int16_t *>(in), qint64(frameCount));
  reinterpret_cast<Application *>(userData)
      ->_decoder->write(static_cast<const int16_t *>(in), qint64(frameCount));
  return paContinue;
}
