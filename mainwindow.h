#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "ui_mainwindow.h"

//// include ctk
//#include <ctkDICOMDatabase.h>
//#include <ctkDICOMIndexer.h>
//#include <ctkFileDialog.h>
//
//include QT
#include <QHash>
#include <QString>
#include <QStringList>
#include <QVariant>
#include <QWidget>
#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class Ui_MainWindow;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
