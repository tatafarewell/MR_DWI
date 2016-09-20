#include "DiffusionCore.h"
#include "ui_DiffusionCore.h"

#include <DicomHelper.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkDICOMreader.h"

// Qt
#include <QCheckBox>
#include <QMessageBox>

DiffusionCore::DiffusionCore(QWidget *parent)
	: QWidget(parent)
	, m_Controls(nullptr)
{
	Initialize();
	CreateQtPartControl(this);
}

DiffusionCore::~DiffusionCore()
{

}

void DiffusionCore::Initialize(){
	//OnImageFilesLoaded(fileLists)
}


void DiffusionCore::CreateQtPartControl(QWidget *parent)
{
	this->m_Controls = new Ui::DiffusionModule;
	this->m_Controls->setupUi(this);

	//vtkSmartPointer< vtkDICOMImageReader > reader =
	//	vtkSmartPointer< vtkDICOMImageReader >::New();
	//reader->SetDirectoryName(argv[1]);
	//reader->Update();
	//int imageDims[3];
	//reader->GetOutput()->GetDimensions(imageDims);


	//for (int i = 0; i < 3; i++)
	//{
	//	riw[i] = vtkSmartPointer< vtkResliceImageViewer >::New();
	//}

	//this->ui->view1->SetRenderWindow(riw[0]->GetRenderWindow());
	//riw[0]->SetupInteractor(
	//	this->ui->view1->GetRenderWindow()->GetInteractor());

	//this->ui->view2->SetRenderWindow(riw[1]->GetRenderWindow());
	//riw[1]->SetupInteractor(
	//	this->ui->view2->GetRenderWindow()->GetInteractor());

	//this->ui->view3->SetRenderWindow(riw[2]->GetRenderWindow());
	//riw[2]->SetupInteractor(
	//	this->ui->view3->GetRenderWindow()->GetInteractor());

}

void DiffusionCore::OnImageFilesLoaded(const QStringList& fileLists)
{
	qDebug() << "Diffusion Module is Loading Files :" << fileLists;
	//vtkStringArray *loadingFiles = fileLists;
}