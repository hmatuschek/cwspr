#include "mainwindow.hh"
#include <QKeyEvent>
#include <QDebug>
#include <QVBoxLayout>
#include <QToolBar>
#include <QComboBox>
#include <QPushButton>
#include <QIntValidator>
#include <QAction>
#include <QLabel>
#include <QMenu>
#include <QActionGroup>
#include <QToolButton>
#include <QSettings>
#include <QInputDialog>
#include <QStatusBar>
#include <QDateTime>


MainWindow::MainWindow(Application &app, QWidget *parent)
    : QMainWindow(parent), _app(app)
{
  setWindowTitle(tr("CWspr by DM3MAT"));
  setMinimumWidth(680);

  QSettings settings;
  ModeId modeId = ModeId(settings.value("mode", MODE_WSPR).toUInt());
  double freq = settings.value("freq", 600.).toDouble();
  QString message = settings.value("message", "DM3MAT JO63 37").toString();
  QString call = settings.value("call", "DM3MAT").toString();
  QString loc  = settings.value("locator", "JO62").toString();
  int dBm      = settings.value("power", 37).toInt();

  QToolBar *toolbar = new QToolBar();
  toolbar->setMovable(true);
  toolbar->setToolButtonStyle(Qt::ToolButtonIconOnly);

  _start = toolbar->addAction(QIcon::fromTheme("media-playback-start", QIcon("://play.png")),
                              tr("Start/Stop"), this, SLOT(onStartToggled(bool)));
  _start->setCheckable(true);
  _start->setChecked(false);
  QMenu *menu = new QMenu();

  menu->addSection(tr("Mode"));
  QMenu *mode = menu->addMenu(tr("Mode"));
  _modes = new QActionGroup(menu);
  _modes->setExclusive(true);
  connect(_modes, SIGNAL(triggered(QAction*)), this, SLOT(onSetMode(QAction*)));

  QAction *wspr = _modes->addAction(tr("WSPR"));
  wspr->setData(uint(MODE_WSPR));
  wspr->setCheckable(true);
  if (MODE_WSPR == modeId)
    wspr->setChecked(true);
  mode->addAction(wspr);
  QAction *jt4 = _modes->addAction(tr("JT4"));
  jt4->setData(uint(MODE_JT4));
  jt4->setCheckable(true);
  if (MODE_JT4 == modeId)
      jt4->setChecked(true);
  mode->addAction(jt4);

  _freq = menu->addAction(tr("Frequency: %1 Hz").arg(int(freq)), this, SLOT(onSetFreq()));
  menu->addSection(tr("Station"));
  _call = menu->addAction(tr("Call: %1").arg(call), this, SLOT(onSetCall()));
  _loc  = menu->addAction(tr("Locator: %1").arg(loc), this, SLOT(onSetLocator()));
  _dBm  = menu->addAction(tr("Power: %1 dBm").arg(dBm), this, SLOT(onSetPower()));

  _menuButton = new QToolButton();
  _menuButton->setText(tr("Settings"));
  _menuButton->setIcon(QIcon::fromTheme("preferences-system", QIcon("://cog.png")));
  _menuButton->setPopupMode(QToolButton::InstantPopup);
  _menuButton->setMenu(menu);
  toolbar->addWidget(_menuButton);

  _message = new QLineEdit(message);
  connect(_message, SIGNAL(returnPressed()), this, SLOT(onSetText()));
  toolbar->addWidget(_message);
  this->addToolBar(Qt::TopToolBarArea, toolbar);

  _plot = new Waterfall(app);
  _plot->setMinimumHeight(150);

  _rx = new QListWidget();

  _ptt = new QLabel();
  _ptt->setText("RX");


  _status = new QStatusBar();
  if (MODE_WSPR == modeId)
    _status->showMessage("WSPR");
  else if (MODE_JT4 == modeId)
    _status->showMessage("JT4");
  _status->addPermanentWidget(_ptt);
  this->setStatusBar(_status);

  QVBoxLayout *layout = new QVBoxLayout();
  layout->setMargin(0); layout->setSpacing(0);
  layout->addWidget(_plot);
  layout->addWidget(_rx);
  QWidget *panel = new QWidget();
  panel->setLayout(layout);

  setCentralWidget(panel);

  connect(&_app, SIGNAL(ptt(bool)), this, SLOT(onPTT(bool)));
  connect(&_app, SIGNAL(newMessage(QString,float,float)), this, SLOT(onNewMessage(QString,float,float)));
}


void
MainWindow::onStartToggled(bool start) {
  if (start) {
    _start->setIcon(QIcon::fromTheme("media-playback-stop", QIcon("://stop.png")));
    /// @todo Verify text for selected encoder.
    _message->setEnabled(false);
    _menuButton->setEnabled(false);
    _app.start();
  } else {
    _start->setIcon(QIcon::fromTheme("media-playback-start", QIcon("://play.png")));
    _app.stop();
    _message->setEnabled(true);
    _menuButton->setEnabled(true);
  }
}

void MainWindow::onSetFreq() {
  QSettings settings;
  double freq = settings.value("freq", 600.).toDouble();
  bool ok=false;
  freq = QInputDialog::getInt(
        nullptr, tr("Set tone-frequency."), tr("CW tone frequency:"),
        int(freq), 300, 2500, 1, &ok);
  if (ok)
    settings.setValue("freq", freq);
  _freq->setText(tr("Frequency: %1 Hz").arg(int(freq)));
}

void
MainWindow::onSetMode(QAction *action) {
  ModeId mode = ModeId(action->data().toUInt());
  QSettings settings;
  settings.setValue("mode", uint(mode));
}

void MainWindow::onSetCall() {
  QSettings settings;
  QString call = settings.value("call", "DM3MAT").toString();
  bool ok=false;
  call = QInputDialog::getText(
        nullptr, tr("Set Call."), tr("Your call:"),
        QLineEdit::Normal, call, &ok);
  if (ok)
    settings.setValue("call", call);
  _call->setText(tr("Call: %1").arg(call));
}

void MainWindow::onSetLocator() {
  QSettings settings;
  QString loc = settings.value("locator", "JO62").toString();
  bool ok=false;
  loc = QInputDialog::getText(
        nullptr, tr("Set Locator."), tr("Your locator:"),
        QLineEdit::Normal, loc, &ok);
  if (ok)
    settings.setValue("locator", loc);
  _loc->setText(tr("Locator: %1").arg(loc));
}

void MainWindow::onSetPower() {
  QSettings settings;
  int pwr = settings.value("power", 37).toInt();
  bool ok=false;
  pwr = QInputDialog::getInt(
        nullptr, tr("Set power."), tr("Transmit power:"),
        int(pwr), 0, 60, 1, &ok);
  if (ok)
    settings.setValue("power", pwr);
  _dBm->setText(tr("Power: %1 dBm").arg(pwr));
}

void
MainWindow::onSetText() {
  QSettings settings;
  settings.setValue("message", _message->text());
}

void
MainWindow::onPTT(bool tx) {
  if (tx)
    _ptt->setText("TX");
  else
    _ptt->setText("RX");
}

void
MainWindow::onNewMessage(QString msg, float freq, float snr) {
  QDateTime now = QDateTime::currentDateTimeUtc();
  msg = QString("%1 @ %2Hz (%3dB): %4").arg(now.toString()).arg(freq).arg(snr).arg(msg);
  _rx->insertItem(0, msg);
}

void
MainWindow::keyPressEvent(QKeyEvent *event) {
  if (Qt::Key_Down == event->key()) {
    _app.key(true);
  }
  if (Qt::Key_Left == event->key()) {
    _app.ditKey(true);
  }
  if (Qt::Key_Right == event->key()) {
    _app.daKey(true);
  }
}


void
MainWindow::keyReleaseEvent(QKeyEvent *event) {
  if (Qt::Key_Down == event->key()) {
    _app.key(false);
  }
  if (Qt::Key_Left == event->key()) {
    _app.ditKey(false);
  }
  if (Qt::Key_Right == event->key()) {
    _app.daKey(false);
  }
}

