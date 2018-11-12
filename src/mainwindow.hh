#ifndef MAINWINDOW_HH
#define MAINWINDOW_HH

#include <QMainWindow>
#include <QLineEdit>
#include <QListView>
#include <QAction>
#include "application.hh"
#include "waterfall.hh"


class MainWindow: public QMainWindow
{
	Q_OBJECT

public:
  MainWindow(Application &app, QWidget *parent=nullptr);

protected slots:
  void onStartToggled(bool start);
  void onSetFreq();
  void onSetMode(QAction *action);
  void onSetCall();
  void onSetLocator();
  void onSetPower();
  void onSetText();
  void onPTT(bool tx);

protected:
  virtual void keyPressEvent(QKeyEvent *event);
  virtual void keyReleaseEvent(QKeyEvent *event);

protected:
  Application &_app;
  QAction *_start;
  QToolButton *_menuButton;
  QActionGroup *_modes;
  QAction *_freq;
  QAction *_call;
  QAction *_loc;
  QAction *_dBm;
  QLineEdit *_message;
  Waterfall *_plot;
  QListView *_rx;
  QStatusBar *_status;
  QLabel *_ptt;
};

#endif // MAINWINDOW_HH
