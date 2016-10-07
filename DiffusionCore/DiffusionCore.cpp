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
#include "QVTKWidget.h"

#include <DicomHelper.h>

// Qt
#include <QCheckBox>
#include <QMessageBox>
#include <QDialog>

//jiangli add
#include <itkVtkAllHeaders.h>

#include <QMouseEvent>

DiffusionCore::DiffusionCore(QWidget *parent)
	:QWidget(parent)
	//, dicomHelp(nullptr)
{	
	this->m_Controls = nullptr;
	this->sourceImage = vtkSmartPointer<vtkImageData>::New();
	this->m_RenderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	m_SourceImageCurrentSlice = 0;
	m_QuantitativeImageCurrentSlice = 0;
	CreateQtPartControl(this);

	this->m_DicomHelper = nullptr;
	this->m_MaskThreshold = 20;
	this->m_ComputedBValue = 2000;
}

DiffusionCore::~DiffusionCore()
{
	//delete m_Controls;
	//this->m_Controls = NULL;

	this->sourceImage->Delete();
	this->sourceImage = NULL;
}

void DiffusionCore::CreateQtPartControl(QWidget *parent)
{
	if (!m_Controls)
	{
		this->m_Controls = new Ui::DiffusionModule;
		this->m_Controls->setupUi(parent);
		m_Controls->bSlider->setMaximum(5000); //maximum supported b value;
		m_Controls->bSlider->setDecimals(0);
		m_Controls->bSlider->setSingleStep(100);
		m_Controls->bSlider->setTickInterval(100);
		m_Controls->bSlider->setValue(2000);
		m_Controls->ThreshSlider->setMaximum(200); //maximum threshhold value;
		m_Controls->ThreshSlider->setDecimals(0);
		m_Controls->ThreshSlider->setSingleStep(5);
		m_Controls->ThreshSlider->setTickInterval(5);
		m_Controls->ThreshSlider->setValue(20);

		//connect Buttons and handle visibility
		
		connect(m_Controls->adcToggle, SIGNAL(toggled(bool)), this, SLOT(calcADC(bool)));		
		connect(m_Controls->pushButton, SIGNAL(toggled(bool)), this, SLOT(calcEADC(bool)));
		connect(m_Controls->cdwiToggle, SIGNAL(toggled(bool)), this, SLOT(calcCDWI(bool)));

		connect(m_Controls->ThreshSlider, SIGNAL(valueChanged(double)), this, SLOT(adc(double)));
		connect(m_Controls->bSlider, SIGNAL(valueChanged(double)), this, SLOT(cDWI(double)));

		//if (m_Controls->adcToggle->isChecked() || m_Controls->eadcButton->isChecked() || m_Controls->cdwiToggle->isChecked())
		//{
		//	qDebug() << "At leaset one of three buttons is checked" << endl;
		//	//connect(m_Controls->eadcButton, SIGNAL(toggled(bool)), this, SLOT(calcEADC(bool)));
		//	//connect(m_Controls->cdwiToggle, SIGNAL(toggled(bool)), this, SLOT(calcCDWI(bool)));
		//	//connect(m_Controls->ThreshLabel->setVisible(true);
		//	m_Controls->ThreshSlider->setVisible(true);
		//	if (m_Controls->cdwiToggle->isChecked())
		//	{
		//		m_Controls->bLabel->setVisible(true);
		//		m_Controls->bSlider->setVisible(true);
		//	}
		//}

		m_Controls->groupBox->hide();
		m_Controls->procBox->hide(); 
		m_Controls->cDWI->hide();
	}

}

void DiffusionCore::OnImageFilesLoaded(const QStringList& fileLists)
{
	
	vtkStringArray* loadingFiles = vtkStringArray::New();
	loadingFiles->SetNumberOfValues(fileLists.size());
	for (int i = 0; i < fileLists.size(); i++)
	{
		loadingFiles->SetValue(i, fileLists.at(i).toStdString().c_str());
	}

	this->m_DicomHelper = new DicomHelper(loadingFiles);
	this->sourceImage = m_DicomHelper->DicomReader->GetOutput();

	DisplayDicomInfo(sourceImage);

	//Tmp used here, add to dicomhelper 
	//int dim[3];
	//this->sourceImage->GetDimensions(dim);
	////m_SourceImageMaxSlice = dim[2] - 1;
	//qDebug() << "B value ="<<dicomHelp.numberOfBValue << endl;

	if (m_DicomHelper->numberOfGradDirection < 2)
	{
		this->m_Controls->faButton->setDisabled(true);
		this->m_Controls->colorFAButton->setDisabled(true);
		this->m_Controls->ivimButton->setDisabled(true);
		this->m_Controls->dtiNameTag->setText("Data does not contain multiple direction, view Only");

		if (m_DicomHelper->numberOfBValue < 2)
		{
			QMessageBox::StandardButton reply;
			reply = QMessageBox::information(this, tr("QMessageBox::information()"), tr("Must Select Diffusion Weighted Image Data"));
			//neutralize all DWI related functions
			this->m_Controls->adcToggle->setDisabled(true);
			this->m_Controls->pushButton->setDisabled(true);
			this->m_Controls->cdwiToggle->setDisabled(true);
			//this->m_Controls->faButton->setDisabled(true);
			//this->m_Controls->colorFAButton->setDisabled(true);
			//this->m_Controls->ivimButton->setDisabled(true);

			this->m_Controls->dataName->setText("Not DWI data, view Only");
			this->m_Controls->dtiNameTag->setText("Not DWI data, view Only");

		}
		else{
			this->m_Controls->dtiNameTag->setText("DTI data loaded");
		}
	}
	else{		
		this->m_Controls->dataName->setText("DWI data loaded");
	}

	QVTKWidget *vtkWindow1 = new QVTKWidget;
	this->m_Controls->displayLayout_2->addWidget(vtkWindow1, 0, 0, 1, 1);	
	SourceImageViewer2D(this->sourceImage, vtkWindow1);
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
	qDebug() << "extent: " << extent[0] << "x" << extent[1] << "x" << extent[2] << "x" << extent[3] << "x" << extent[4] << "x" << extent[5] << endl;
	imageData->GetScalarRange(range);
	qDebug() << "range: " << range[0] << "x" << range[1] << endl;
	
	std::cout << " scaleSlope = " << this->m_DicomHelper->scaleSlope << std::endl;
	std::cout << "scaleIntercept = " << this->m_DicomHelper->scaleIntercept << std::endl;
	std::cout << "----------------------diffusion related parameters--------------------------------" << std::endl;

	for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i = i + this->m_DicomHelper->numberOfGradDirection)
	{
		std::cout << "bValueList " << i << ": " << this->m_DicomHelper->BvalueList.at(i / this->m_DicomHelper->numberOfGradDirection) << std::endl;
	}
	std::cout << "numberOfComponents = " << this->m_DicomHelper->numberOfComponents << std::endl;
	//std::cout << "numberOfComponents = " << this->m_DicomHelper->numberOfComponents << std::endl;
	std::cout << "numberOfGradDirection = " << this->m_DicomHelper->numberOfGradDirection << std::endl;
	std::cout << "numberOfBValue = " << this->m_DicomHelper->numberOfBValue << std::endl;
}

void DiffusionCore::SourceImageViewer2D(vtkSmartPointer <vtkImageData> imageData, QVTKWidget *qvtkWidget)
{
	double *imageDataRange = new double[2];
	imageDataRange = imageData->GetScalarRange();//Replace with to be displayed

	vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
	imageViewer->SetInputData(imageData);
	imageViewer->SetSliceOrientationToXY();
	imageViewer->SetColorWindow(imageDataRange[1] - imageDataRange[0]);
	imageViewer->SetColorLevel(0.5* (imageDataRange[1] + imageDataRange[0]));

	vtkSmartPointer<vtkTextProperty> sliceTextProp = vtkSmartPointer<vtkTextProperty>::New();
	sliceTextProp->SetFontFamilyToCourier();
	sliceTextProp->SetFontSize(18);
	sliceTextProp->SetVerticalJustificationToBottom();
	sliceTextProp->SetJustificationToLeft();

	vtkSmartPointer<vtkTextMapper> sliceTextMapper = vtkSmartPointer<vtkTextMapper>::New();
	std::string msg = StatusMessage::Format(imageViewer->GetSliceMin(), imageViewer->GetSliceMax());
	sliceTextMapper->SetInput(msg.c_str());
	sliceTextMapper->SetTextProperty(sliceTextProp);

	vtkSmartPointer<vtkActor2D> sliceTextActor = vtkSmartPointer<vtkActor2D>::New();
	sliceTextActor->SetMapper(sliceTextMapper);
	sliceTextActor->SetPosition(15, 10);

	imageViewer->GetRenderer()->AddActor2D(sliceTextActor);

	//vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle = vtkSmartPointer<myVtkInteractorStyleImage>::New();
	//vtkSmartPointer<vtkTestCallbackCommand> testCallbackCommand = vtkSmartPointer<vtkTestCallbackCommand>::New();

	myInteractorStyle->SetImageViewer(imageViewer);
	myInteractorStyle->SetStatusMapper(sliceTextMapper);
	myInteractorStyle->GetCurrentSliceNumber(m_SourceImageCurrentSlice);
	//imageViewer->SetupInteractor(renderWindowInteractor);
	//renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	//renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, testCallbackCommand);
	imageViewer->SetupInteractor(m_RenderWindowInteractor);
	m_RenderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	imageViewer->GetRenderWindow()->SetSize(qvtkWidget->width(), qvtkWidget->height());
	//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
	qvtkWidget->SetRenderWindow(imageViewer->GetRenderWindow());
	//connect(qvtkWidget, SIGNAL(mouseEvent(QMouseEvent*);), this, SLOT(outputMOUSEevent(QMouseEvent*)));
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetSize(200, 200);
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetPosition(this->x(), this->y());
	//qvtkWidget->GetRenderWindow()->GetInteractor()->SetSize(200, 200);

	qvtkWidget->show();
	std::cout << " qvtkWidget height " << qvtkWidget->height() << " qvtkWidget width " << qvtkWidget->width() << std::endl;
	std::cout << " Now in SourceImage viewer " << std::endl;
	imageViewer->GetRenderer()->ResetCamera(); //Reset camera and then render is better
	imageViewer->Render();
	//renderWindowInteractor->Initialize();
	//renderWindowInteractor->Start();
	m_RenderWindowInteractor->Initialize();
	m_RenderWindowInteractor->Start();

	// Needed?
	delete[] imageDataRange;
	imageDataRange = NULL;
	//deleter qvtkWidget somewhere!!! @Yang Ming
	//imageData is SmartPointer, it should know when to delete
}

void DiffusionCore::QuantitativeImageViewer2D(vtkSmartPointer <vtkImageData> imageData, QVTKWidget *qvtkWidget)
{
	double *imageDataRange = new double[2];
	imageDataRange = imageData->GetScalarRange();	

	vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
	imageViewer->SetInputData(imageData);
	imageViewer->SetSliceOrientationToXY();
	imageViewer->SetColorWindow(imageDataRange[1] - imageDataRange[0]);
	imageViewer->SetColorLevel(0.5* (imageDataRange[1] + imageDataRange[0]));

	vtkSmartPointer<vtkTextProperty> sliceTextProp = vtkSmartPointer<vtkTextProperty>::New();
	sliceTextProp->SetFontFamilyToCourier();
	sliceTextProp->SetFontSize(18);
	sliceTextProp->SetVerticalJustificationToBottom();
	sliceTextProp->SetJustificationToLeft();

	vtkSmartPointer<vtkTextMapper> sliceTextMapper = vtkSmartPointer<vtkTextMapper>::New();
	std::string msg = StatusMessage::Format(imageViewer->GetSliceMin(), imageViewer->GetSliceMax());
	sliceTextMapper->SetInput(msg.c_str());
	sliceTextMapper->SetTextProperty(sliceTextProp);

	vtkSmartPointer<vtkActor2D> sliceTextActor = vtkSmartPointer<vtkActor2D>::New();
	sliceTextActor->SetMapper(sliceTextMapper);
	sliceTextActor->SetPosition(15, 10);

	imageViewer->GetRenderer()->AddActor2D(sliceTextActor);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle = vtkSmartPointer<myVtkInteractorStyleImage>::New();

	myInteractorStyle->SetImageViewer(imageViewer);
	myInteractorStyle->SetStatusMapper(sliceTextMapper);
	myInteractorStyle->GetCurrentSliceNumber(m_QuantitativeImageCurrentSlice);
	imageViewer->SetupInteractor(renderWindowInteractor);
	renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	imageViewer->GetRenderWindow()->SetSize(qvtkWidget->width(),qvtkWidget->height());
	//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
	qvtkWidget->SetRenderWindow(imageViewer->GetRenderWindow());
	qvtkWidget->show();

	std::cout << " Now in Quantitative ImageViewer " << std::endl;
	imageViewer->GetRenderer()->ResetCamera(); //Reset camera and then render is better
	imageViewer->Render();
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

	delete[] imageDataRange;
	imageDataRange = NULL;
}

void DiffusionCore::calcADC(bool _istoggled)
{	
	if (_istoggled)
	{

		
		//Set Threshhold bar visible
		this->m_Controls->groupBox->setVisible(_istoggled);
		//Visualize data
		QVTKWidget *vtkWindow2 = new QVTKWidget;
		this->m_Controls->displayLayout_2->addWidget(vtkWindow2, 0, 1, 1, 1);


		//===================================================================================================
		//Will Pack below into a filter
		typedef itk::Image < unsigned short, 3> SrcImageType;
		typedef itk::Image < float, 3> FloatImageType;
		typedef itk::VectorContainer< unsigned int, FloatImageType::Pointer > ImageContainerType;

		typedef itk::VTKImageToImageFilter <SrcImageType>	VtkToItkConverterType;
		typedef itk::CastImageFilter< SrcImageType, FloatImageType >	CastFilterType;
		typedef itk::ShiftScaleImageFilter<FloatImageType, FloatImageType>	ShiftScaleType;


		ImageContainerType::Pointer imageContainer = ImageContainerType::New();
		imageContainer->Reserve(this->m_DicomHelper->numberOfComponents);
		std::cout << "m_DicomHelper numof components" << this->m_DicomHelper->numberOfComponents << std::endl;
		for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i++)
		{
			//Handle each scalar component indivisually
			vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
			scalarComponent->SetInputData(this->sourceImage);
			scalarComponent->SetComponents(i);
			scalarComponent->Update();//Crutial, otherwise abort after running

			//VTK to ITK Image Data
			VtkToItkConverterType::Pointer vtkToItkImageFilter = VtkToItkConverterType::New();
			vtkToItkImageFilter->SetInput(scalarComponent->GetOutput());
			//vtkToItkImageFilter->Update();

			//unsigned short image to float image
			CastFilterType::Pointer castFilter = CastFilterType::New();
			castFilter->SetInput(vtkToItkImageFilter->GetOutput());
			//castFilter->Update();

			//Shift and scale signal back to FP value
			//Take some time to finish the computation
			ShiftScaleType::Pointer shiftScale = ShiftScaleType::New();
			shiftScale->SetInput(castFilter->GetOutput());
			shiftScale->SetShift(-m_DicomHelper->scaleIntercept);
			shiftScale->SetScale(1.0 / m_DicomHelper->scaleSlope);
			shiftScale->Update();

			//Save vector image & image container
			imageContainer->InsertElement(i, dynamic_cast <FloatImageType*> (shiftScale->GetOutput()));
		}

		//Get vector Image
		typedef itk::ComposeImageFilter<FloatImageType>		ImageToVectorImageType;
		ImageToVectorImageType::Pointer imageToVectorImageFilter = ImageToVectorImageType::New();

		for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i++)
		{
			imageToVectorImageFilter->SetInput(i, imageContainer->GetElement(i));
		}

		imageToVectorImageFilter->Update();
		//std::cout << "vectorImage: vectorLength = " << imageToVectorImageFilter->GetOutput()->GetVectorLength() << std::endl;

		//////////////////////////
		// Sort Data order
		//B0 image should be first!!!
		//data order: b0, b11, b12...b13.... bij....bn1, bn2, ...bnm (n = bvalues, m = directions)


		//Mask image filter
		typedef itk::MaskVectorImageFilter <float> MaskFilterType;
		MaskFilterType::Pointer maskFilter = MaskFilterType::New();
		maskFilter->SetInput(imageToVectorImageFilter->GetOutput());//input is a Image Pointer!!!
		maskFilter->SetMaskThreshold(m_MaskThreshold / m_DicomHelper->scaleSlope);//Get from UI or user-interaction
		maskFilter->Update();//output is a Image Pointer!!!
		//std::cout << "maskFilter: vectorLength = " << maskFilter->GetOutput()->GetVectorLength() << std::endl;

		std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
		// Adc Map filter
		typedef itk::AdcMapFilter <float, float> AdcMapFilterType;
		AdcMapFilterType::Pointer adcMap = AdcMapFilterType::New();
		adcMap->SetInput(maskFilter->GetOutput());
		adcMap->SetBValueList(this->m_DicomHelper->BvalueList);
		adcMap->Update();
		std::cout << "------------- AdcMapFilter end runing ----------- " << std::endl;
		//std::cout << "adcVectorImage: vectorLength = " << adcMap->GetOutput()->GetVectorLength() << std::endl;

		//cDwi filter
		typedef itk::ComputedDwiFilter <float, float> ComputedDwiFilterType;
		ComputedDwiFilterType::Pointer computedDwi = ComputedDwiFilterType::New();
		computedDwi->SetInput(adcMap->GetOutput());
		computedDwi->SetNumOfDiffDirections(this->m_DicomHelper->numberOfGradDirection);
		computedDwi->SetComputedBValue(m_ComputedBValue);//Get from UI input
		computedDwi->Update();
		//std::cout << "cDWi vectorImage: vectorLength = " << computedDwi->GetOutput()->GetVectorLength() << std::endl;

		//vector image to scalar image or imageContainer
		typedef itk::VectorIndexSelectionCastImageFilter <itk::VectorImage<float, 3>, FloatImageType> VectorImageToImageType;
		VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
		vectorImageToImageFilter->SetIndex(0);
		vectorImageToImageFilter->SetInput(adcMap->GetOutput());
		vectorImageToImageFilter->Update();


		//Data Clipping
		typedef itk::DisplayOptimizer < FloatImageType, FloatImageType> DisplayOptimizerType;
		DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
		displayOptimizer->SetInput(vectorImageToImageFilter->GetOutput());
		displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
		displayOptimizer->Update();

		//Rescale signal intensity to display
		typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
		RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
		rescaleFilter->SetInput(displayOptimizer->GetOutput());
		rescaleFilter->SetOutputMaximum(2000.0);// Lucky Guess? 
		rescaleFilter->SetOutputMinimum(0.0);
		rescaleFilter->Update();
		std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		///////////////////////////////////////////
		//ITK to VTK for visualization
		typedef itk::ImageToVTKImageFilter <FloatImageType> itkToVtkConverter;
		itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
		convItkToVtk->SetInput(rescaleFilter->GetOutput());
		convItkToVtk->Update();

		
		//Will Pack above into a filter
		//===================================================================================================

		this->cacheImage = convItkToVtk->GetOutput();
		//debug out

		QuantitativeImageViewer2D(this->cacheImage, vtkWindow2);
	}
	else{
		//Hide Threshhold bar
		this->m_Controls->groupBox->setVisible(_istoggled);
		//clear up the widget rendering ADC		
		QLayoutItem *existItem = this->m_Controls->displayLayout_2->itemAtPosition(0,1); 
		QWidget * existWidget = existItem->widget(); 
		if (existWidget != NULL) {
			this->m_Controls->displayLayout_2->removeWidget(existWidget);
			existWidget->setParent(NULL);//if you want to delete the widget, do: widget->setParent(NULL); delete widget; 
			delete existWidget;
		}
		this->m_Controls->displayLayout_2->update();
	}
	
}
 
void DiffusionCore::calcEADC(bool _istoggled)
{

}

void DiffusionCore::calcCDWI(bool _istoggled)
{
	if (_istoggled)
	{
		this->m_Controls->cDWI->setVisible(_istoggled);
		//Visualize data
		QVTKWidget *vtkWindow3 = new QVTKWidget;
		this->m_Controls->displayLayout_2->addWidget(vtkWindow3, 1, 0, 1, 1);

		
			//===================================================================================================
			//Will Pack below into a filter
			typedef itk::Image < unsigned short, 3> SrcImageType;
			typedef itk::Image < float, 3> FloatImageType;
			typedef itk::VectorContainer< unsigned int, FloatImageType::Pointer > ImageContainerType;

			typedef itk::VTKImageToImageFilter <SrcImageType>	VtkToItkConverterType;
			typedef itk::CastImageFilter< SrcImageType, FloatImageType >	CastFilterType;
			typedef itk::ShiftScaleImageFilter<FloatImageType, FloatImageType>	ShiftScaleType;


			ImageContainerType::Pointer imageContainer = ImageContainerType::New();
			imageContainer->Reserve(this->m_DicomHelper->numberOfComponents);
			std::cout << "m_DicomHelper numof components" << this->m_DicomHelper->numberOfComponents << std::endl;
	vtkSmartPointer <vtkExtractVOI> ExtractVOI = vtkSmartPointer <vtkExtractVOI>::New();
	ExtractVOI->SetInputData(this->sourceImage);
	ExtractVOI->SetVOI(0, this->m_DicomHelper->imageDimensions[0] - 1, 0, this->m_DicomHelper->imageDimensions[1] - 1, m_SourceImageCurrentSlice, m_SourceImageCurrentSlice);
	ExtractVOI->Update();
	std::cout << "---------------------------- VOI is correct ?--------------------" << std::endl;
	DisplayDicomInfo(ExtractVOI->GetOutput());
			for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i++)
			{
				//Handle each scalar component indivisually
				vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
		scalarComponent->SetInputConnection(ExtractVOI->GetOutputPort());
				scalarComponent->SetComponents(i);
				scalarComponent->Update();//Crutial, otherwise abort after running

				//VTK to ITK Image Data
				VtkToItkConverterType::Pointer vtkToItkImageFilter = VtkToItkConverterType::New();
				vtkToItkImageFilter->SetInput(scalarComponent->GetOutput());
				//vtkToItkImageFilter->Update();

				//unsigned short image to float image
				CastFilterType::Pointer castFilter = CastFilterType::New();
				castFilter->SetInput(vtkToItkImageFilter->GetOutput());
				//castFilter->Update();

				//Shift and scale signal back to FP value
				//Take some time to finish the computation
				ShiftScaleType::Pointer shiftScale = ShiftScaleType::New();
				shiftScale->SetInput(castFilter->GetOutput());
				shiftScale->SetShift(-m_DicomHelper->scaleIntercept);
				shiftScale->SetScale(1.0 / m_DicomHelper->scaleSlope);
				shiftScale->Update();

				//Save vector image & image container
				imageContainer->InsertElement(i, dynamic_cast <FloatImageType*> (shiftScale->GetOutput()));
			}

			//Get vector Image
			typedef itk::ComposeImageFilter<FloatImageType>		ImageToVectorImageType;
			ImageToVectorImageType::Pointer imageToVectorImageFilter = ImageToVectorImageType::New();

			for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i++)
			{
				imageToVectorImageFilter->SetInput(i, imageContainer->GetElement(i));
			}

			imageToVectorImageFilter->Update();
			//std::cout << "vectorImage: vectorLength = " << imageToVectorImageFilter->GetOutput()->GetVectorLength() << std::endl;

			//////////////////////////
			// Sort Data order
			//B0 image should be first!!!
			//data order: b0, b11, b12...b13.... bij....bn1, bn2, ...bnm (n = bvalues, m = directions)


			//Mask image filter
			typedef itk::MaskVectorImageFilter <float> MaskFilterType;
			MaskFilterType::Pointer maskFilter = MaskFilterType::New();
			maskFilter->SetInput(imageToVectorImageFilter->GetOutput());//input is a Image Pointer!!!
			maskFilter->SetMaskThreshold(m_MaskThreshold / m_DicomHelper->scaleSlope);//Get from UI or user-interaction
			maskFilter->Update();//output is a Image Pointer!!!
			//std::cout << "maskFilter: vectorLength = " << maskFilter->GetOutput()->GetVectorLength() << std::endl;

			std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
			// Adc Map filter
			typedef itk::AdcMapFilter <float, float> AdcMapFilterType;
			AdcMapFilterType::Pointer adcMap = AdcMapFilterType::New();
			adcMap->SetInput(maskFilter->GetOutput());
			adcMap->SetBValueList(this->m_DicomHelper->BvalueList);
			adcMap->Update();
			std::cout << "------------- AdcMapFilter end runing ----------- " << std::endl;
			//std::cout << "adcVectorImage: vectorLength = " << adcMap->GetOutput()->GetVectorLength() << std::endl;

			//cDwi filter
			typedef itk::ComputedDwiFilter <float, float> ComputedDwiFilterType;
			ComputedDwiFilterType::Pointer computedDwi = ComputedDwiFilterType::New();
			computedDwi->SetInput(adcMap->GetOutput());
			computedDwi->SetNumOfDiffDirections(this->m_DicomHelper->numberOfGradDirection);
			computedDwi->SetComputedBValue(m_ComputedBValue);//Get from UI input
			computedDwi->Update();
			//std::cout << "cDWi vectorImage: vectorLength = " << computedDwi->GetOutput()->GetVectorLength() << std::endl;

			//vector image to scalar image or imageContainer
			typedef itk::VectorIndexSelectionCastImageFilter <itk::VectorImage<float, 3>, FloatImageType> VectorImageToImageType;
			VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
			vectorImageToImageFilter->SetIndex(0);
			vectorImageToImageFilter->SetInput(computedDwi->GetOutput());
			vectorImageToImageFilter->Update();


			//Data Clipping
			typedef itk::DisplayOptimizer < FloatImageType, FloatImageType> DisplayOptimizerType;
			DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
			displayOptimizer->SetInput(vectorImageToImageFilter->GetOutput());
			displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
			displayOptimizer->Update();

			//Rescale signal intensity to display
			typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
			RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
			rescaleFilter->SetInput(displayOptimizer->GetOutput());
			rescaleFilter->SetOutputMaximum(2000.0);// Lucky Guess? 
			rescaleFilter->SetOutputMinimum(0.0);
			rescaleFilter->Update();
			std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
			std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

			///////////////////////////////////////////
			//ITK to VTK for visualization
			typedef itk::ImageToVTKImageFilter <FloatImageType> itkToVtkConverter;
			itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
			convItkToVtk->SetInput(rescaleFilter->GetOutput());
			convItkToVtk->Update();
		


		//Will Pack above into a filter
		//===================================================================================================

		QuantitativeImageViewer2D(convItkToVtk->GetOutput(), vtkWindow3);
	}
	else{
		this->m_Controls->bLabel->hide();
		this->m_Controls->bSlider->hide();
		//clear up the widget rendering ADC		
		QLayoutItem *existItem = this->m_Controls->displayLayout_2->itemAtPosition(1, 0);
		QWidget * existWidget = existItem->widget();
		if (existWidget != NULL) {
			this->m_Controls->displayLayout_2->removeWidget(existWidget);
			existWidget->setParent(NULL);//if you want to delete the widget, do: widget->setParent(NULL); delete widget; 
			delete existWidget;
		}
		this->m_Controls->displayLayout_2->update();
	}

}

//void DiffusionCore::ADCmethod()
//{
//}


void DiffusionCore::cDWI(double computedBValue)
{
	//qDebug() << "computedBValue from slider is " << computedBValue << endl;
	m_ComputedBValue = computedBValue;
}

void DiffusionCore::adc(double maskThreshold)
{
	QLayoutItem *existItem = this->m_Controls->displayLayout_2->itemAtPosition(0, 1);
	QVTKWidget * existWidget = static_cast<QVTKWidget* >(existItem->widget());
	m_MaskThreshold = maskThreshold;
	//ImageViewer2D(ADCmethod(), existWidget);
	//qDebug() << "maskThreshold value from slider is " << maskThreshold << endl;

}