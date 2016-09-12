#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "ui_mainhub.h"

//// include ctk
//#include <ctkDICOMDatabase.h>
//#include <ctkDICOMIndexer.h>
//#include <ctkFileDialog.h>

//#include "QmitkDicomExternalDataWidget.h"
//#include "QmitkDicomLocalStorageWidget.h"
//
//include QT
#include <QHash>
#include <QString>
#include <QStringList>
#include <QVariant>
#include <QWidget>
#include <QMainWindow>

class Ui_MainWindow;

class MainHub : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainHub(QWidget *parent = 0);
    ~MainHub();

private:
	Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
