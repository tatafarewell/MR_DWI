#include <QApplication>
#include "mainhub.h"

int main( int argc, char** argv )
{
  // QT Stuff
    QApplication app( argc, argv );

    MainHub myUI;
	myUI.show();

    return app.exec();
}
