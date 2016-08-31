#include <QApplication>
#include "mainwindow.h"

int main( int argc, char** argv )
{
  // QT Stuff
    QApplication app( argc, argv );

    MainWindow myUI;
	myUI.show();

    return app.exec();
}
