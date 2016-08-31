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
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}
