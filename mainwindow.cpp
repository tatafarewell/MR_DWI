#include "mainwindow.h"
#include "ui_mainwindow.h"

//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include "vtkResliceImageViewer.h"
//#include "vtkResliceCursorLineRepresentation.h"
//#include "vtkResliceCursorThickLineRepresentation.h"
//#include "vtkResliceCursorWidget.h"
//#include "vtkResliceCursorActor.h"
//#include "vtkResliceCursorPolyDataAlgorithm.h"
//#include "vtkResliceCursor.h"
//#include "vtkDICOMImageReader.h"
//#include "vtkCellPicker.h"
//#include "vtkProperty.h"
//#include "vtkPlane.h"
//#include "vtkImageData.h"
//#include "vtkCommand.h"
//#include "vtkPlaneSource.h"
//#include "vtkLookupTable.h"
//#include "vtkImageMapToWindowLevelColors.h"
//#include "vtkInteractorStyleImage.h"
//#include "vtkImageSlabReslice.h"
//#include "vtkBoundedPlanePointPlacer.h"
//#include "vtkDistanceWidget.h"
//#include "vtkDistanceRepresentation.h"
//#include "vtkHandleRepresentation.h"
//#include "vtkResliceImageViewerMeasurements.h"
//#include "vtkDistanceRepresentation2D.h"
//#include "vtkPointHandleRepresentation3D.h"
//#include "vtkPointHandleRepresentation2D.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
	ui(new Ui::MainWindow)
{
	//Initialize();
	//CreateQtPartControl(this);
	ui->setupUi(this);
	connect(ui->externalDataWidget, SIGNAL(SignalStartDicomImport(QStringList)), ui->diffusionModule, SLOT(OnImageFilesLoaded(QStringList)));
}

//void MainWindow::CreateQtPartControl(QWidget *parent)
//{
//	// build up qt Widget, unless already done
//	if (!ui)
//	{
//		ui = new Ui::MainWindow;
//		ui->setupUi(this);
//		ui->viewExternalDataButton->setVisible(true);
//		ui->ctkDICOMtable->setTableOrientation(Qt::Vertical);
//		ui->ctkDICOMtable->setDICOMDatabase(m_ExternalDatabase);
//
//		this->SetupImportDialog();
//		this->SetupProgressDialog(parent);
//
//		//connect Buttons
//		connect(ui->downloadButton, SIGNAL(clicked()), this, SLOT(OnDownloadButtonClicked()));
//		connect(ui->viewExternalDataButton, SIGNAL(clicked()), this, SLOT(OnViewButtonClicked()));
//		connect(ui->directoryButton, SIGNAL(clicked()), this, SLOT(OnScanDirectory()));
//
//		connect(ui->ctkDICOMtable, SIGNAL(seriesSelectionChanged(const QStringList&)),this, SLOT(OnSeriesSelectionChanged(const QStringList&)));
//		connect(ui->ctkDICOMtable, SIGNAL(seriesDoubleClicked(const QModelIndex&)),this, SLOT(OnViewButtonClicked()));
//
//		connect(m_ProgressDialog, SIGNAL(canceled()), m_ExternalIndexer, SLOT(cancel()));
//		connect(m_ExternalIndexer, SIGNAL(indexingComplete()), this, SLOT(OnFinishedImport()));
//
//		connect(m_ExternalIndexer, SIGNAL(indexingFilePath(const QString&)), m_ProgressDialogLabel, SLOT(setText(const QString&)));
//		connect(m_ExternalIndexer, SIGNAL(progress(int)), m_ProgressDialog, SLOT(setValue(int)));
//	}
//}

//void MainWindow::Initialize()
//{
//	//m_ExternalDatabase = new ctkDICOMDatabase(this);
//
//	////try{
//	//	m_ExternalDatabase->openDatabase(QString(":memory:"), QString("EXTERNAL-DB"));
//	////}
//	////catch (std::exception e){
//	////	std::cout<< "Database Open Error";
//	////	m_ExternalDatabase->closeDatabase();
//	////	return;
//	////}
//
//	//m_ExternalIndexer = new ctkDICOMIndexer(this);
//
//
//}

MainWindow::~MainWindow()
{
    delete ui;
}
