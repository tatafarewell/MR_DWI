#include "DiffusionCore.h"
#include "ui_DiffusionCore.h"

#include <vtkDataObjectToTable.h>
#include <vtkElevationFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkQtTableView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkVectorText.h>
#include <vtkStringArray.h>
#include <vtkDICOMReader.h>
#include <vtkImageData.h>
#include "vtkImageViewer2.h"
#include "vtkSmartPointer.h"

#include <DicomHelper.h>

// Qt
#include <QCheckBox>
#include <QMessageBox>
#include <QDialog>

DiffusionCore::DiffusionCore(QWidget *parent)
	: QWidget(parent)
	, m_Controls(nullptr)
	, sourceImage(vtkSmartPointer<vtkImageData>::New())
	//, dicomHelp(nullptr)
{
	CreateQtPartControl(this);
}

DiffusionCore::~DiffusionCore()
{

}

//void DiffusionCore::Initialize(){
//	//OnImageFilesLoaded(fileLists)
//}
//
//
void DiffusionCore::CreateQtPartControl(QWidget *parent)
{
	if (!m_Controls)
	{
		this->m_Controls = new Ui::DiffusionModule;
		this->m_Controls->setupUi(parent);
		m_Controls->bSlider->setMaximum(5000); //maximum supported b value;
		m_Controls->ThreshSlider->setMaximum(200); //maximum threshhold value;
		//connect Buttons
		connect(m_Controls->pushButton, SIGNAL(clicked()), this, SLOT(calcADC()));
		connect(m_Controls->bSlider, SIGNAL(valueChanged(int)), this, SLOT(cDWI(int)));
		connect(m_Controls->ThreshSlider, SIGNAL(valueChanged(int)), this, SLOT(adc(int)));
	}
}

void  DiffusionCore::calcADC(){}

void DiffusionCore::cDWI(int bvalue)
{
}

void DiffusionCore::adc(int threshhold)
{
}


void DiffusionCore::OnImageFilesLoaded(const QStringList& fileLists)
{
	
	vtkStringArray* loadingFiles = vtkStringArray::New();
	loadingFiles->SetNumberOfValues(fileLists.size());
	for (int i = 0; i < fileLists.size(); i++)
	{
		loadingFiles->SetValue(i, fileLists.at(i).toStdString().c_str());
	}

	//vtkSmartPointer<vtkDICOMReader> DicomReader = vtkSmartPointer<vtkDICOMReader>::New();
	//DicomReader->SetFileNames(loadingFiles);
	//DicomReader->Update();

	DicomHelper dicomHelp(loadingFiles);
	sourceImage = dicomHelp.DicomReader->GetOutput();

	DisplayDicomInfo(sourceImage);
	qDebug() << "B value ="<<dicomHelp.numberOfBValue << endl;

	if (dicomHelp.numberOfBValue < 1)
	{

		QMessageBox::StandardButton reply;
		reply = QMessageBox::information(this, tr("QMessageBox::information()"), tr("Must Select Diffusion Weighted Image Data"));

		//QMessageBox msgBox(QMessageBox::Warning, tr("QMessageBox::warning()"),
		//tr("Must Select Diffusion Weighted Image Data"), 0, this);
		//msgBox.setDetailedText(tr("Diffusion Weighted Data has to have b value > 0"));
		//msgBox.addButton(tr("Save &Again"), QMessageBox::AcceptRole);
		//msgBox.addButton(tr("&Understand"), QMessageBox::RejectRole); 
		//if (msgBox.exec() == QMessageBox::AcceptRole)			
		//else
	}

	vtkSmartPointer<vtkImageViewer2> viewer1 = vtkSmartPointer<vtkImageViewer2>::New();

	this->m_Controls->qvtkWidget->SetRenderWindow( viewer1->GetRenderWindow());
	viewer1->SetupInteractor( this->m_Controls->qvtkWidget->GetRenderWindow()->GetInteractor());
	viewer1->SetInputData(sourceImage);
	viewer1->SetSliceOrientation(0);
	this->m_Controls->qvtkWidget->show();
	//viewer1->SetResliceModeToAxisAligned();

}

void DiffusionCore::DisplayDicomInfo(vtkSmartPointer <vtkImageData> imageData)
{

	//cout << "fileNames: " << reader->GetFileNames()<< endl;
	qDebug() << "number of COMPONENTS: " << imageData->GetNumberOfScalarComponents() << endl;
	const int dataDim = imageData->GetDataDimension();
	//int cells = reader->GetOutput()->GetNumberOfCells();
	//int points = reader->GetOutput()->GetNumberOfPoints();
	//cout << "points: " << points << endl;
	//cout << "cells: " << cells << endl;
	int dims[3];
	double origins[3];
	double spacing[3];
	int extent[6];
//	double center[3];
	double range[2];
	imageData->GetDimensions(dims);
	qDebug() << "image dims: " << dims[0] << "x" << dims[1] << "x" << dims[2] << endl;
	imageData->GetOrigin(origins);
	qDebug() << "image origins: " << origins[0] << " " << origins[1] << " " << origins[2] << endl;
	imageData->GetSpacing(spacing);
	qDebug() << "image spacing: " << spacing[0] << "x" << spacing[1] << "x" << spacing[2] << endl;
	imageData->GetExtent(extent);
	qDebug() << "extent: " << extent[0] << "x" << extent[1] << "x" << extent[2] << endl;
	imageData->GetScalarRange(range);
	qDebug() << "range: " << range[0] << "x" << range[1] << endl;
}