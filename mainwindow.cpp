#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "DicomModule.h"

#include <QMessageBox>
#include <QStyleFactory>
#include <qfile.h>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
	ui(new Ui::MainWindow)
{
	//Initialize();
	//CreateQtPartControl(this);
	DicomUI = new DicomModule(this);
	DicomUI->hide();

	//Set Window Style
	QFile file("style.qss");
	file.open(QFile::ReadOnly);
	QString styleSheet = QLatin1String(file.readAll());
	setStyleSheet(styleSheet);

	ui->setupUi(this);
	showMaximized();

	//ui->toolLine->setPalette(framePalette);
	//ui->toolLine->setAutoFillBackground(true); // set dockwidget as trasparent floating
	//ui->dockWidget->setWindowFlags(Qt::FramelessWindowHint);
	//ui->dockWidget->setFloating(true);

	QApplication::setStyle(QStyleFactory::create("fusion"));
	connect(ui->FileButton, SIGNAL(clicked()), this, SLOT(onStartdicom()));
	connect(DicomUI, SIGNAL(SignalDicomRead(QStringList)), ui->diffusionModule, SLOT(OnImageFilesLoaded(QStringList)));
	connect(ui->diffusionModule, SIGNAL(SignalDicomLoaded(bool)), this, SLOT(onDicomLoaded(bool)));
	//connect(DicomUI, SIGNAL(SignalDicomToDataManager(QStringList)), ui->diffusionModule, SLOT(OnImageFilesLoaded(QStringList)));
	//connect(ui->DicomUI, SIGNAL(SignalStartDicomImport(QStringList)), ui->internalDataWidget, SLOT(OnStartDicomImport(QStringList)));
}

void MainWindow::onStartdicom()
{
	//std::cout << "button clicked" << std::endl;
	DicomUI->setWindowFlags(Qt::Window);
	DicomUI->setWindowTitle(tr("DICOM BROWSER"));
	DicomUI->resize(1600, 1000);
	DicomUI->show();
    //connect();
}

void MainWindow::onDicomLoaded(bool isLoaded)
{
	//std::cout << "received SignalDicomLoaded" << std::endl;
	if (isLoaded)
	{
		DicomUI->hide();
	}
	
}

MainWindow::~MainWindow()
{
    delete ui;
}
