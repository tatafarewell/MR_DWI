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

//Customized Header
#include <itkVtkAllHeaders.h>
#include <itkGetDiffusionImageFilter.h>
#include <itkTensor.h>
#include "vtkImageChangeInformation.h"
#include <VTKtoITK.h>
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
		m_Controls->bSlider->setTracking(false);
		m_Controls->ThreshSlider->setMaximum(200); //maximum threshhold value;
		m_Controls->ThreshSlider->setDecimals(0);
		m_Controls->ThreshSlider->setSingleStep(5);
		m_Controls->ThreshSlider->setTickInterval(5);
		m_Controls->ThreshSlider->setValue(20);
		m_Controls->ThreshSlider->setTracking(false);
		//connect Buttons and handle visibility

		this->m_Controls->ADCTool->setDisabled(true);
		this->m_Controls->DTITool->setDisabled(true);

		//this->m_Controls->adcToggle->setDisabled(true);
		//this->m_Controls->eadcToggle->setDisabled(true);
		//this->m_Controls->cdwiToggle->setDisabled(true);
		//this->m_Controls->faToggle->setDisabled(true);
		//this->m_Controls->colorFAToggle->setDisabled(true);
		//this->m_Controls->ivimToggle->setDisabled(true);

		connect(m_Controls->adcToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcADC(bool)));
		connect(m_Controls->eadcToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcEADC(bool)));
		connect(m_Controls->cdwiToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcCDWI(bool)));
		connect(m_Controls->faToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcFA(bool)));
		connect(m_Controls->colorFAToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcColorFA(bool)));
		connect(m_Controls->ivimToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcIVIM(bool)));

		connect(m_Controls->ThreshSlider, SIGNAL(valueChanged(double)), this, SLOT(onThreshSlide(double)));
		//connect(m_Controls->ThreshSlider, SIGNAL(valueChanged(double)), this, SLOT(onCalcADC(true)));//jiangli test
		connect(m_Controls->bSlider, SIGNAL(valueChanged(double)), this, SLOT(onBSlide(double)));
		//m_Controls->DTITool->hide();
		//m_Controls->ADCTool->hide();
		m_Controls->cDWI->hide();
		m_Controls->Thresh->hide();
	}
}

void DiffusionCore::OnImageFilesLoaded(const QStringList& fileLists)
{	
	//enable all the buttons 
	this->m_Controls->ADCTool->setEnabled(true);
	this->m_Controls->DTITool->setEnabled(true);

	//this->m_Controls->adcToggle->setEnabled(true);
	//this->m_Controls->eadcToggle->setEnabled(true);
	//this->m_Controls->cdwiToggle->setEnabled(true);
	//this->m_Controls->faToggle->setEnabled(true);
	//this->m_Controls->colorFAToggle->setEnabled(true);
	//this->m_Controls->ivimToggle->setEnabled(true);

	vtkStringArray* loadingFiles = vtkStringArray::New();
	loadingFiles->SetNumberOfValues(fileLists.size());
	for (int i = 0; i < fileLists.size(); i++)
	{
		loadingFiles->SetValue(i, fileLists.at(i).toStdString().c_str());
	}

	this->m_DicomHelper = new DicomHelper(loadingFiles);
	this->sourceImage = m_DicomHelper->DicomReader->GetOutput();

	//DisplayDicomInfo(sourceImage);


	//Tmp used here, add to dicomhelper 
	//int dim[3];
	//this->sourceImage->GetDimensions(dim);
	////m_SourceImageMaxSlice = dim[2] - 1;
	//qDebug() << "B value ="<<dicomHelp.numberOfBValue << endl;

	if (!m_DicomHelper->tensorComputationPossible)
	{
		this->m_Controls->faToggle->setDisabled(true);
		this->m_Controls->colorFAToggle->setDisabled(true);		
		//this->m_Controls->dtiNameTag->setText("Data does not contain multiple direction, view Only");

		if (m_DicomHelper->numberOfBValue < 2)
		{
			QMessageBox::StandardButton reply;
			reply = QMessageBox::information(this, tr("QMessageBox::information()"), tr("Must Select Diffusion Weighted Image Data"));
			//neutralize all DWI related functions
			this->m_Controls->adcToggle->setDisabled(true);
			this->m_Controls->eadcToggle->setDisabled(true);
			this->m_Controls->cdwiToggle->setDisabled(true);
			if (m_DicomHelper->numberOfBValue < 3)
			{
				this->m_Controls->ivimToggle->setDisabled(true);
			}
			//this->m_Controls->faToggle->setDisabled(true);
			//this->m_Controls->colorFAToggle->setDisabled(true);
			//this->m_Controls->ivimToggle->setDisabled(true);
		}
		else{			
		}
	}
	else{		

	}
	
	QVTKWidget *vtkWindow1 = new QVTKWidget;
	this->m_Controls->displayLayout->addWidget(vtkWindow1, 0, 0, 1, 1);	
	//update the layout table
	std::vector<int> row1; row1.push_back(0);
	layoutTable.push_back(row1);

	std::cout << "Emitting signal" << std::endl;
	emit SignalDicomLoaded(true);

	SourceImageViewer2D(this->sourceImage, vtkWindow1);

}

void DiffusionCore::onCalcADC(bool _istoggled) //It's a SLOT!!!
{

	if (_istoggled)
	{
		//Set Threshhold bar visible
		this->m_Controls->Thresh->setVisible(_istoggled);
		
		QVTKWidget *vtkWindow2 = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow2, ADC);

		//this->m_Controls->displayLayout->addWidget(vtkWindow2, 0, 1, 1, 1);
		
		//layoutTable() = 1;

		vtkSmartPointer <vtkImageData> calculatedAdc = vtkSmartPointer <vtkImageData>::New();
		this->AdcCalculator(calculatedAdc);

		QuantitativeImageViewer2D(calculatedAdc, vtkWindow2,"ADC");

	}
	else
	{
		//Hide Threshhold bar
		this->m_Controls->Thresh->setVisible(_istoggled);

		//initialize first: clear up the widget 
		//Or should we just clear the renderer in current window?
		ui_RemoveWindow(ADC);

	}
	
}

void DiffusionCore::AdcCalculator(vtkSmartPointer <vtkImageData> imageData)
{
	/////////////////////////////////////
	//preprocessing: from VTK to ITK vector images
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
	std::cout << "---------------------------- VOI is correct ? ---------------------" << std::endl;
	this->DisplayDicomInfo(ExtractVOI->GetOutput());


	for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i++)
	{
		//Handle each scalar component indivisually
		vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
		scalarComponent->SetInputData(ExtractVOI->GetOutput());
		scalarComponent->SetComponents(i);
		scalarComponent->Update();//Crutial, otherwise abort after running

		//VTK to ITK Image Data
		VtkToItkConverterType::Pointer vtkToItkImageFilter = VtkToItkConverterType::New();
		vtkToItkImageFilter->SetInput(scalarComponent->GetOutput());
		vtkToItkImageFilter->Update();

		//unsigned short image to float image
		CastFilterType::Pointer castFilter = CastFilterType::New();
		castFilter->SetInput(vtkToItkImageFilter->GetOutput());
		castFilter->Update();

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

	//vector image to scalar image or imageContainer
	typedef itk::VectorIndexSelectionCastImageFilter <itk::VectorImage<float, 3>, FloatImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
	vectorImageToImageFilter->SetIndex(1);
	vectorImageToImageFilter->SetInput(adcMap->GetOutput());
	vectorImageToImageFilter->Update();


	//Rescale signal intensity to display
	typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
	rescaleFilter->SetInput(vectorImageToImageFilter->GetOutput());
	rescaleFilter->SetOutputMaximum(2000.0);
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->Update();
	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	//Data Clipping
	typedef itk::DisplayOptimizer < FloatImageType, SrcImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
	displayOptimizer->SetInput(rescaleFilter->GetOutput());
	displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer->Update();

	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	///////////////////////////////////////////
	//ITK to VTK
	typedef itk::ImageToVTKImageFilter <SrcImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
	convItkToVtk->SetInput(displayOptimizer->GetOutput());
	convItkToVtk->Update();

	imageData->DeepCopy(convItkToVtk->GetOutput());
	//Crutial not imageData = convItkToVtk->GetOutput(),
	//Because convItkToVtk->Update() pointer is recycled after this caller!!!

	//Visualize data
	//QVTKWidget *vtkComputedDwiWindow = new QVTKWidget;
	//this->m_Controls->displayLayout->addWidget(vtkComputedDwiWindow, 0, 1, 1, 1);
	//this->m_Controls->displayLayout->update();
	//QuantitativeImageViewer2D(convItkToVtk->GetOutput(), vtkComputedDwiWindow);
}

void DiffusionCore::onCalcEADC(bool _istoggled) //It's a SLOT!!!
{
	if (_istoggled)
	{
		this->m_Controls->Thresh->setVisible(_istoggled);

		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, EADC);

		vtkSmartPointer <vtkImageData> calculatedEAdc = vtkSmartPointer <vtkImageData>::New();
		this->EAdcCalculator(calculatedEAdc);
		QuantitativeImageViewer2D(calculatedEAdc, vtkWindow,"eADC");

	}
	else
	{
		this->m_Controls->Thresh->setVisible(_istoggled);
		
		//clear up existing widget
		ui_RemoveWindow(EADC);
	}
}

void DiffusionCore::onCalcCDWI(bool _istoggled)  
{


	if (_istoggled)
	{
		this->m_Controls->cDWI->setVisible(_istoggled);		

		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, CDWI);
		
		vtkSmartPointer <vtkImageData> computedDwi = vtkSmartPointer <vtkImageData>::New();
		this->ComputedDwi(computedDwi);
		QuantitativeImageViewer2D(computedDwi, vtkWindow,"cDWI");
	}
	else
	{
		this->m_Controls->cDWI->setVisible(_istoggled);		
		//clear up existing widget
		ui_RemoveWindow(CDWI);
	}
}

void DiffusionCore::ComputedDwi(vtkSmartPointer <vtkImageData> imageData)
{
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


	//Rescale signal intensity to display
	typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
	rescaleFilter->SetInput(vectorImageToImageFilter->GetOutput());
	rescaleFilter->SetOutputMaximum(2000.0);
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->Update();
	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	//Data Clipping
	typedef itk::DisplayOptimizer < FloatImageType, SrcImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
	displayOptimizer->SetInput(rescaleFilter->GetOutput());
	displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer->Update();

	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	///////////////////////////////////////////
	//ITK to VTK for visualization
	typedef itk::ImageToVTKImageFilter < SrcImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
	convItkToVtk->SetInput(displayOptimizer->GetOutput());
	convItkToVtk->Update();

	imageData->DeepCopy(convItkToVtk->GetOutput());
}

void DiffusionCore::onThreshSlide(double maskThreshold) //It's a SLOT!!!
{
	m_MaskThreshold = maskThreshold;
	DiffusionCore::ShareWindowEvent();
}

void DiffusionCore::onBSlide(double computedBValue) //It's a SLOT!!!
{
	m_ComputedBValue = computedBValue;

	//clear up existing widget
	QLayoutItem *existItem = this->m_Controls->displayLayout->itemAtPosition(1, 0);
	QWidget * existWidget = existItem->widget();
	if (existWidget != NULL) 
	{
		QVTKWidget *cDwiWidget= static_cast <QVTKWidget*> (existWidget);
		vtkSmartPointer <vtkImageData> computedDwi = vtkSmartPointer <vtkImageData>::New();
		this->ComputedDwi(computedDwi);
		QuantitativeImageViewer2D(computedDwi, cDwiWidget, "cDWI");
	}
}

void DiffusionCore::onCalcFA(bool _istoggled)
{
	if (_istoggled)
	{
		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, FA);
		
		vtkSmartPointer <vtkImageData> calculatedFA = vtkSmartPointer <vtkImageData>::New();
		this->FaCalculator(calculatedFA);

		QuantitativeImageViewer2D(calculatedFA, vtkWindow, "FA");

	}
	else
	{
		ui_RemoveWindow(FA);
	}
}

void DiffusionCore::onCalcColorFA(bool _istoggled)
{
	if (_istoggled)
	{
		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, CFA);

		vtkSmartPointer <vtkImageData> calculatedColorFA = vtkSmartPointer <vtkImageData>::New();
		this->ColorFACalculator(calculatedColorFA);
		QuantitativeImageViewer2D(calculatedColorFA, vtkWindow, "Color FA" );
	}
	else
	{
		ui_RemoveWindow(CFA);
	}
}

void DiffusionCore::onCalcIVIM(bool _istoggled)
{
	if (_istoggled)
	{
		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, IVIM);
	}
	else
	{
		ui_RemoveWindow(IVIM);
	}
}

void DiffusionCore::ui_InsertWindow(int& rowInd, int& colInd, QVTKWidget *vtkWindow, imageType imageLabel)
{
	//int i(0),winInd(0);
	bool flag(true);
	int cur_RowNum = layoutTable.size();
	int cur_ColNum = layoutTable[0].size();
	std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>";
	std::cout << "cur girdlayout is " << cur_RowNum << " Row " << cur_ColNum << " Col " << std::endl;

	//Default case: a rect shape is filled, fill the the new window at either the first pos of a new row or a new col
	if (cur_RowNum < cur_ColNum)
	{
		rowInd = cur_RowNum;
		colInd = 0;
		//update layoutTable accordingly, add one row
		std::vector<int> aNewRow(cur_ColNum, NOIMAGE);
	}
	else{
		rowInd = 0;
		colInd = cur_ColNum;
	}

	//Extra case: a rect shape is not filled, Note that it always fills column first. hence only loop the columns of the last row.
	for (int i = 0; i < cur_ColNum&&flag; i++)
	{
		QLayoutItem *existItem = m_Controls->displayLayout->itemAtPosition(cur_RowNum - 1, i);
		if (!existItem)
		{
			rowInd = cur_RowNum - 1 ;//last Row
			colInd = i;
			flag = false;
			std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>it is not fully filled" << std::endl;
		}
	}
	
	std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>";
	std::cout << "Will add "<<imageLabel<<" window to [ " << rowInd << ":" << colInd << "]" << std::endl;

	this->m_Controls->displayLayout->addWidget(vtkWindow, rowInd, colInd, 1, 1);	


	//Update layOutTable
	if (colInd < 1 && rowInd >0 ) //if a new line is added
	{
		std::vector<int> aNewRow(1, imageLabel);
		layoutTable.push_back(aNewRow);
	}
	else{
		layoutTable[rowInd].push_back(imageLabel);
	}

	//debug info starts
	std::vector<std::vector<int>>::iterator rowit;
	std::vector<int>::iterator colit;
	for (rowit = layoutTable.begin(); rowit != layoutTable.end(); rowit++)
	{
		std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>";
		for (colit = (*rowit).begin(); colit != (*rowit).end(); ++colit)
		{
			std::cout << *colit << " ";
		}
		std::cout << std::endl;
	}
	//debug info ends

}

void DiffusionCore::ui_RemoveWindow(imageType imageLabel)
{
	int colWidth = layoutTable[0].size();
	int del_rowInd(0), del_colInd(0);	
	int rowCounter(0), colCounter(0);
	int cur_RowNum = layoutTable.size();
	int cur_ColNum = layoutTable[0].size();

	std::vector<std::vector<int>>::iterator rowit;
	std::vector<int>::iterator colit;

	//STEP1, Find window rendering assigned image and delete it, update the layoutTable. 
	for (rowit = layoutTable.begin(),rowCounter=0; rowit != layoutTable.end(); rowit++, rowCounter++)
	{
		for (colit = (*rowit).begin(),colCounter=0; colit != (*rowit).end();)
		{
			//std::cout << "Running at " << rowCounter << ":" << colCounter << std::endl;
			if (*colit == imageLabel)
			{
				colit = (*rowit).erase(colit);
				del_rowInd = rowCounter;
				del_colInd = colCounter;
			}
			else
			{
				++colit; colCounter++;
			}
		}
	}	

	std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>deleting the widget type "<< imageLabel <<" at "<<del_rowInd<<":"<<del_colInd << std:: endl;
	
	for (rowit = layoutTable.begin(); rowit != layoutTable.end(); rowit++)
	{
		std::cout << "[GRIDLAYOUT DEBUG INFO] BEFORE REORDER WINDOW>>>>>>>> |";
		for (colit = (*rowit).begin(); colit != (*rowit).end(); ++colit)
		{
			std::cout << *colit << " ";
		}
		std::cout << "|" << std::endl;
	}

	//STEP2, Re-order the layoutTable and delete correct window
	if (ui_IsWdWSquare())
	{
		std::cout << "[GRIDLAYOUT DEBUG INFO] >>>>>>>> SQUARE REORDER" << cur_ColNum << std::endl;
		ui_dumpWindow(0, cur_ColNum-1);//delete the last window in first row
		std::cout << "[GRIDLAYOUT DEBUG INFO] >>>>>>>> WDW DUMPED" << std::endl;
		if (del_rowInd != 0)//If deleting element not at first row.
		{
			layoutTable[del_rowInd].push_back(*(layoutTable[0].end() - 1));
			layoutTable[0].erase(layoutTable[0].end() - 1);
		}
		
	}
	else{
		std::cout << "[GRIDLAYOUT DEBUG INFO] >>>>>>>> NON-SQUARE REORDER" << std::endl;

		if (del_rowInd == cur_RowNum - 1)
		{			
			ui_dumpWindow(cur_RowNum - 1, layoutTable[cur_RowNum - 1].size());
			std::cout << "[GRIDLAYOUT DEBUG INFO] >>>>>>>> WDW DUMPED" << std::endl;
			if (layoutTable[del_rowInd].size() == 0)
			{
				layoutTable.erase(layoutTable.end() - 1); //delete the last row if it has no element
			}

		}
		else
		{
			layoutTable[del_rowInd].push_back(*(layoutTable[cur_RowNum-1].end() - 1)); //put the last element to the del row
			layoutTable[cur_RowNum-1].erase(layoutTable[cur_RowNum-1].end() - 1);
			ui_dumpWindow(cur_RowNum - 1, layoutTable[cur_RowNum - 1].size());
			std::cout << "[GRIDLAYOUT DEBUG INFO] >>>>>>>> WDW DUMPED" << std::endl;
			if (layoutTable[cur_RowNum - 1].size() == 0)
			{
				layoutTable.erase(layoutTable.end() - 1); //delete the last row if it has no element
			}
		}

	}

	for (rowit = layoutTable.begin(); rowit != layoutTable.end(); rowit++)
	{
		std::cout << "[GRIDLAYOUT DEBUG INFO] AFTER REORDER WINDOW>>>>>>>> |";
		for (colit = (*rowit).begin(); colit != (*rowit).end(); ++colit)
		{
			std::cout << *colit << " ";
		}
		std::cout <<"|"<< std::endl;
	}

	//STEP3, Then according to layoutTable to re-render all images. 
	
	ui_reDrawAllWdw();
	
	std::cout << "[GRIDLAYOUT DEBUG INFO] After Rendering girdlayout is ";
	std::cout << m_Controls->displayLayout->rowCount() << " Row " << m_Controls->displayLayout->columnCount() << " Col " << std::endl;

	this->m_Controls->displayLayout->update();	

}

bool DiffusionCore::ui_IsWdWSquare()
{
	float WdWNum(0.0);
	std::vector<std::vector<int>>::iterator rowit;
	for (rowit = layoutTable.begin(); rowit != layoutTable.end(); rowit++)
	{
		WdWNum += (*rowit).size();
	}

	if (sqrt(WdWNum) == floor(sqrt(WdWNum)))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void DiffusionCore::ui_dumpWindow(int row, int col)
{
	QLayoutItem *existItem = this->m_Controls->displayLayout->itemAtPosition(row, col);
	QWidget * existWidget = existItem->widget();
	if (existWidget != NULL) {
		std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>deleting window " << " at " << row << ":" << col << std::endl;
		this->m_Controls->displayLayout->removeWidget(existWidget);
		existWidget->setParent(NULL);//if you want to delete the widget, do: widget->setParent(NULL); delete widget; 
		delete existWidget;
	}
	this->m_Controls->displayLayout->update();
	//std::cout << "[GRIDLAYOUT DEBUG INFO] After Deleting girdlayout is ";
	//std::cout << m_Controls->displayLayout->rowCount() << " Row " << m_Controls->displayLayout->columnCount() << " Col " << std::endl;
}

void DiffusionCore::ui_reDrawAllWdw()
{
	std::cout << "[GRIDLAYOUT DEBUG INFO] Before Rendering girdlayout is ";
	std::cout << layoutTable.size() << " X " << layoutTable[0].size() << std::endl;

	for (int rowCounter = 0; rowCounter < layoutTable.size(); rowCounter++)
	{
		for (int colCounter = 0; colCounter < layoutTable[rowCounter].size(); colCounter++)
		{
			std::cout << rowCounter << "-" << layoutTable[rowCounter].size() << "-" << (colCounter != 0 || rowCounter != 0) << std::endl;
			if (colCounter != 0 || rowCounter != 0)
			{
				std::cout << "Rendering Image Type" << layoutTable[rowCounter][colCounter] << "at Layout" << rowCounter << "-" << colCounter << std::endl;
				QLayoutItem *existItem = this->m_Controls->displayLayout->itemAtPosition(rowCounter, colCounter);
				QWidget * existWidget = existItem->widget();
				if (existWidget != NULL)
				{
					QVTKWidget *thisWindow = static_cast <QVTKWidget*> (existWidget);
					vtkSmartPointer <vtkImageData> thisImage = vtkSmartPointer <vtkImageData>::New();

					switch (layoutTable[rowCounter][colCounter])
					{
					case ORIGINAL:
						break;
					case ADC:
						this->ComputedDwi(thisImage);
						QuantitativeImageViewer2D(thisImage, thisWindow, "ADC");
						break;
					case CDWI:
						this->AdcCalculator(thisImage);
						QuantitativeImageViewer2D(thisImage, thisWindow, "cDWI");
						break;
					case EADC:
						this->EAdcCalculator(thisImage);
						QuantitativeImageViewer2D(thisImage, thisWindow, "eADC");
						break;
					case FA:
						this->FaCalculator(thisImage);
						QuantitativeImageViewer2D(thisImage, thisWindow, "FA");
						break;
					case CFA:
						this->ColorFACalculator(thisImage);
						QuantitativeImageViewer2D(thisImage, thisWindow, "color FA");
						break;
					case IVIM:
						std::cout << "[GRIDLAYOUT DEBUG INFO] >> should render IVIM" << std::endl;
						//this->IVIMCalculator(thisImage);
						//QuantitativeImageViewer2D(thisImage, thisWindow);
						break;
					case NOIMAGE:
						std::cout << "[GRIDLAYOUT DEBUG INFO] >> NO image was rendered HERE" << std::endl;
						break;
					}
				}
			}
		}
	}
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
	imageData->GetScalarRange(range);//1. cannot be type of float here, it's a bug of vtk?  2. error while calculating quantitative output images: pipeline should be updated before calling this method!!!
	qDebug() << "range: " << range[0] << "x" << range[1] << endl;
	std::cout << " imageData: GetScalarType() = " << imageData->GetScalarType() << std::endl;
	std::cout << " scaleSlope = " << this->m_DicomHelper->scaleSlope << std::endl;
	std::cout << "scaleIntercept = " << this->m_DicomHelper->scaleIntercept << std::endl;

	std::cout << "diffusion related parameters---" << std::endl;
	for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i = i + this->m_DicomHelper->numberOfGradDirection)
	{
		int j = 0;
		std::cout << "bValueList " << j << ": " << this->m_DicomHelper->BvalueList.at(i / this->m_DicomHelper->numberOfGradDirection) << std::endl;
		j++;
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

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//Use QVTKInteractor rather than vtkRenderWindowInteractor!!! So that interactor start and end events are handled by qApp->exec() and qApp->Exit();
	//vtkSmartPointer<QVTKInteractor> renderWindowInteractor = vtkSmartPointer<QVTKInteractor>::New();
	vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle = vtkSmartPointer<myVtkInteractorStyleImage>::New();
	//vtkSmartPointer<vtkTestCallbackCommand> testCallbackCommand = vtkSmartPointer<vtkTestCallbackCommand>::New();

	myInteractorStyle->SetImageViewer(imageViewer);
	myInteractorStyle->SetStatusMapper(sliceTextMapper);
	myInteractorStyle->GetCurrentSliceNumber(m_SourceImageCurrentSlice);
	//imageViewer->SetupInteractor(renderWindowInteractor);
	renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	renderWindowInteractor->AddObserver(vtkCommand::MouseWheelForwardEvent, this, &DiffusionCore::ShareWindowEvent);
	renderWindowInteractor->AddObserver(vtkCommand::MouseWheelBackwardEvent, this, &DiffusionCore::ShareWindowEvent);
	//imageViewer->GetRenderWindow()->SetSize(qvtkWidget->width(), qvtkWidget->height());
	//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
	qvtkWidget->SetRenderWindow(imageViewer->GetRenderWindow());
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetSize(qvtkWidget->width(), qvtkWidget->height());
	qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetSize(800, 800);
	qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetPosition(qvtkWidget->x(), qvtkWidget->y());
	qvtkWidget->GetRenderWindow()->SetInteractor(renderWindowInteractor);//crutial to let qvtkWidget share the same interactor with imageViewer
	qvtkWidget->show();
	//std::cout << " qvtkWidget height " << qvtkWidget->height() << " qvtkWidget width " << qvtkWidget->width() << std::endl;
	std::cout << " Now in SourceImage viewer " << std::endl;
	imageViewer->GetRenderer()->ResetCamera(); //Reset camera and then render is better
	imageViewer->Render();
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();
}

void DiffusionCore::ShareWindowEvent()
{
	ui_reDrawAllWdw();
}

void DiffusionCore::QuantitativeImageViewer2D(vtkSmartPointer <vtkImageData> imageData, QVTKWidget *qvtkWidget, std::string imageLabel)
{
	if (qvtkWidget->GetRenderWindow()->GetInteractor())
	{
		qvtkWidget->GetRenderWindow()->Finalize();
		qvtkWidget->GetRenderWindow()->GetInteractor()->ExitCallback();
	}
	double *imageDataRange = new double[2];
	imageDataRange = imageData->GetScalarRange();//Replace with to be displayed
	
	double colorWindow, colorLevel;
	if (imageData->GetNumberOfScalarComponents() == 3)
	{
		//color map 
		colorWindow = 255.0;
		colorLevel = 127.5;

	}
	else
	{
		double *imageDataRange = new double[2];
		imageDataRange = imageData->GetScalarRange();//Replace with to be displayed
		colorWindow = imageDataRange[1] - imageDataRange[2];
		colorLevel = 0.5* (imageDataRange[1] + imageDataRange[0]);
	}

	vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
	imageViewer->SetInputData(imageData);
	imageViewer->SetSliceOrientationToXY();
	imageViewer->SetColorWindow(colorWindow);
	imageViewer->SetColorLevel(colorLevel);

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

	vtkSmartPointer<vtkTextMapper> imageLabelMapper = vtkSmartPointer<vtkTextMapper>::New();
	//std::string msg = StatusMessage::Format(imageViewer->GetSliceMin(), imageViewer->GetSliceMax());
	imageLabelMapper->SetInput(imageLabel.c_str());
	imageLabelMapper->SetTextProperty(sliceTextProp);

	vtkSmartPointer<vtkActor2D>  imageLabelActor = vtkSmartPointer<vtkActor2D>::New();
	imageLabelActor->SetMapper(imageLabelMapper);
	imageLabelActor->SetPosition(0, 10);

	imageViewer->GetRenderer()->AddActor2D(imageLabelActor);

	//vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//Use QVTKInteractor rather than vtkRenderWindowInteractor!!! So that interactor start and end events are handled by qApp->exec() and qApp->Exit();
	vtkSmartPointer<QVTKInteractor> renderWindowInteractor = vtkSmartPointer<QVTKInteractor>::New();
	vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle = vtkSmartPointer<myVtkInteractorStyleImage>::New();

	myInteractorStyle->SetImageViewer(imageViewer);
	myInteractorStyle->SetStatusMapper(sliceTextMapper);
	myInteractorStyle->GetCurrentSliceNumber(m_QuantitativeImageCurrentSlice);

	renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
	//imageViewer->GetRenderWindow()->SetSize(qvtkWidget->width(),qvtkWidget->height());
	//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
	qvtkWidget->SetRenderWindow(imageViewer->GetRenderWindow());
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetSize(qvtkWidget->width(), qvtkWidget->height());
	qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetSize(800, 800);
	qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetPosition(qvtkWidget->x(), qvtkWidget->y());
	qvtkWidget->GetRenderWindow()->SetInteractor(renderWindowInteractor);//crutial to let qvtkWidget share the same interactor with imageViewer
	//imageViewer->SetupInteractor(renderWindowInteractor);
	qvtkWidget->show();	
	std::cout << " Now in Quantitative ImageViewer " << std::endl;
	imageViewer->GetRenderer()->ResetCamera(); //Reset camera and then render is better
	imageViewer->Render();
	renderWindowInteractor->Initialize();

	//std::cout << "image data quantitative viewer!!!!!!! " << imageViewer->GetInput()->GetScalarRange()[0] << imageViewer->GetInput()->GetScalarRange()[1] << std::endl;
	//renderWindowInteractor->Initialize();
	//renderWindowInteractor->Start();
}
void DiffusionCore::FaCalculator(vtkSmartPointer <vtkImageData> imageData)
{
	/////////////////////////////////////
	//preprocessing: from VTK to ITK vector images
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
	//std::cout << "---------------------------- VOI is correct ? ---------------------" << std::endl;
	//this->DisplayDicomInfo(ExtractVOI->GetOutput());
	//std::cout << "---------------------------- End of VOI ---------------------" << std::endl;

	//PR: mask image doesn't change for slices other than the first one whlie sliding the maskthreshold slider
	//fix extent range here, default z extent is m_SourceImageCurrentSlice to m_SourceImageCurrentSlice;
	//Update image extent, otherwise maskFilter won't work for slices other than the first one (because seed starts at origin)
	//Update the origin as well
	vtkSmartPointer <vtkImageChangeInformation> changeInfo = vtkSmartPointer <vtkImageChangeInformation>::New();
	changeInfo->SetInputData(ExtractVOI->GetOutput());
	changeInfo->SetOutputOrigin(0, 0, 0);
	changeInfo->SetExtentTranslation(0, 0, -m_SourceImageCurrentSlice);
	changeInfo->Update();
	//std::cout << "---------------------------- translateExtent is correct ? ---------------------" << std::endl;
	//this->DisplayDicomInfo(changeInfo->GetOutput());
	//std::cout << "---------------------------- End of translateExtent ---------------------" << std::endl;



	if (m_DicomHelper->tensorComputationPossible)
	{
		//convert data from VTK to ITK: transfer to float and scaling to None scaled data.
		vtkSmartPointer <vtkImageData> imageSliceData = vtkSmartPointer <vtkImageData>::New();
		imageSliceData = changeInfo->GetOutput();
		VTKtoITK convertedToITK(m_DicomHelper->IsoImageLabel,
			m_DicomHelper->scaleSlope, m_DicomHelper->scaleIntercept, m_DicomHelper->numberOfComponents);
		convertedToITK.Scaling(imageSliceData);

		//Mask image filter
		typedef itk::MaskVectorImageFilter <float> MaskFilterType;
		MaskFilterType::Pointer maskFilter = MaskFilterType::New();
		maskFilter->SetInput(convertedToITK.outComposedImage->GetOutput());//input is a Image Pointer!!!
		maskFilter->SetMaskThreshold(m_MaskThreshold / m_DicomHelper->scaleSlope);//Get from UI or user-interaction
		maskFilter->Update();//output is a Image Pointer!!!
		//std::cout << "maskFilter: vectorLength = " << maskFilter->GetOutput()->GetVectorLength() << std::endl;

		//DTI claculation
		cout << "DTI is calculating.." << endl;
		typedef itk::GetDiffusionImageFilter<float, float> GetDiffusionImageType;
		GetDiffusionImageType::Pointer GetDiffusionImages = GetDiffusionImageType::New();
		GetDiffusionImages->SetInput(maskFilter->GetOutput());
		GetDiffusionImages->SetHMatrix(m_DicomHelper->finalH);
		GetDiffusionImages->SetSlice2PatMatrix(m_DicomHelper->slice2PatMatrix);
		GetDiffusionImages->SetBValueList(m_DicomHelper->BvalueList);
		GetDiffusionImages->Update();

		//Adc
		typedef itk::VectorIndexSelectionCastImageFilter <itk::VectorImage<float, 3>, FloatImageType> VectorImageToImageType;
		VectorImageToImageType::Pointer vectorImageToImageFilterAdc = VectorImageToImageType::New();
		vectorImageToImageFilterAdc->SetIndex(2);
		vectorImageToImageFilterAdc->SetInput(GetDiffusionImages->GetOutput1());
		vectorImageToImageFilterAdc->Update();

		//Rescale signal intensity to display
		typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
		RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
		rescaleFilter->SetInput(vectorImageToImageFilterAdc->GetOutput());
		rescaleFilter->SetOutputMaximum(4095.0);
		rescaleFilter->SetOutputMinimum(0.0);
		rescaleFilter->Update();
		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		//Data Clipping
		typedef itk::DisplayOptimizer < FloatImageType, SrcImageType> DisplayOptimizerType;
		DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
		displayOptimizer->SetInput(rescaleFilter->GetOutput());
		displayOptimizer->SetCoveragePercent(0.9);//Default is 0.99
		displayOptimizer->Update();

		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		///////////////////////////////////////////
		//ITK to VTK
		typedef itk::ImageToVTKImageFilter <SrcImageType> itkToVtkConverter;
		itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
		convItkToVtk->SetInput(displayOptimizer->GetOutput());
		convItkToVtk->Update();

		imageData->DeepCopy(convItkToVtk->GetOutput());

	}
	else
	{
		for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i++)
		{
			//Handle each scalar component indivisually
			vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
			scalarComponent->SetInputData(changeInfo->GetOutput());
			scalarComponent->SetComponents(i);
			scalarComponent->Update();//Crutial, otherwise abort after running

			//VTK to ITK Image Data
			VtkToItkConverterType::Pointer vtkToItkImageFilter = VtkToItkConverterType::New();
			vtkToItkImageFilter->SetInput(scalarComponent->GetOutput());
			vtkToItkImageFilter->Update();

			//unsigned short image to float image
			CastFilterType::Pointer castFilter = CastFilterType::New();
			castFilter->SetInput(vtkToItkImageFilter->GetOutput());
			castFilter->Update();

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

		//get vector Image
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

		//std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
		// Adc Map filter
		typedef itk::AdcMapFilter <float, float> AdcMapFilterType;
		AdcMapFilterType::Pointer adcMap = AdcMapFilterType::New();
		adcMap->SetInput(maskFilter->GetOutput());
		adcMap->SetBValueList(this->m_DicomHelper->BvalueList);
		adcMap->Update();

		//std::cout << "------------- AdcMapFilter end runing ----------- " << std::endl;
		//std::cout << "adcVectorImage: vectorLength = " << adcMap->GetOutput()->GetVectorLength() << std::endl;

		//vector image to scalar image or imageContainer
		typedef itk::VectorIndexSelectionCastImageFilter <itk::VectorImage<float, 3>, FloatImageType> VectorImageToImageType;
		VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
		vectorImageToImageFilter->SetIndex(1);
		vectorImageToImageFilter->SetInput(adcMap->GetOutput());
		vectorImageToImageFilter->Update();

		//Rescale signal intensity to display
		typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
		RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
		rescaleFilter->SetInput(vectorImageToImageFilter->GetOutput());
		rescaleFilter->SetOutputMaximum(4095.0);
		rescaleFilter->SetOutputMinimum(0.0);
		rescaleFilter->Update();
		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		//Data Clipping
		typedef itk::DisplayOptimizer < FloatImageType, SrcImageType> DisplayOptimizerType;
		DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
		displayOptimizer->SetInput(rescaleFilter->GetOutput());
		displayOptimizer->SetCoveragePercent(0.9);//Default is 0.99
		displayOptimizer->Update();

		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		///////////////////////////////////////////
		//ITK to VTK
		typedef itk::ImageToVTKImageFilter <SrcImageType> itkToVtkConverter;
		itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
		convItkToVtk->SetInput(displayOptimizer->GetOutput());
		convItkToVtk->Update();

		imageData->DeepCopy(convItkToVtk->GetOutput());
		//Crutial not imageData = convItkToVtk->GetOutput(),
		//Because convItkToVtk->Update() pointer is recycled after this caller!!!

		//Visualize data
		//QVTKWidget *vtkComputedDwiWindow = new QVTKWidget;
		//this->m_Controls->displayLayout_2->addWidget(vtkComputedDwiWindow, 0, 1, 1, 1);
		//this->m_Controls->displayLayout_2->update();
		//QuantitativeImageViewer2D(convItkToVtk->GetOutput(), vtkComputedDwiWindow);
	}


}
void DiffusionCore::ColorFACalculator(vtkSmartPointer <vtkImageData> imageData)
{
	
	/////////////////////////////////////
	//preprocessing: from VTK to ITK vector images
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
	//std::cout << "---------------------------- VOI is correct ? ---------------------" << std::endl;
	//this->DisplayDicomInfo(ExtractVOI->GetOutput());
	//std::cout << "---------------------------- End of VOI ---------------------" << std::endl;

	//PR: mask image doesn't change for slices other than the first one whlie sliding the maskthreshold slider
	//fix extent range here, default z extent is m_SourceImageCurrentSlice to m_SourceImageCurrentSlice;
	//Update image extent, otherwise maskFilter won't work for slices other than the first one (because seed starts at origin)
	//Update the origin as well
	vtkSmartPointer <vtkImageChangeInformation> changeInfo = vtkSmartPointer <vtkImageChangeInformation>::New();
	changeInfo->SetInputData(ExtractVOI->GetOutput());
	changeInfo->SetOutputOrigin(0, 0, 0);
	changeInfo->SetExtentTranslation(0, 0, -m_SourceImageCurrentSlice);
	changeInfo->Update();
	//std::cout << "---------------------------- translateExtent is correct ? ---------------------" << std::endl;
	//this->DisplayDicomInfo(changeInfo->GetOutput());
	//std::cout << "---------------------------- End of translateExtent ---------------------" << std::endl;


		//convert data from VTK to ITK: transfer to float and scaling to None scaled data.
		vtkSmartPointer <vtkImageData> imageSliceData = vtkSmartPointer <vtkImageData>::New();
		imageSliceData = changeInfo->GetOutput();
		VTKtoITK convertedToITK(m_DicomHelper->IsoImageLabel,
			m_DicomHelper->scaleSlope, m_DicomHelper->scaleIntercept, m_DicomHelper->numberOfComponents);
		convertedToITK.Scaling(imageSliceData);

		//Mask image filter
		typedef itk::MaskVectorImageFilter <float> MaskFilterType;
		MaskFilterType::Pointer maskFilter = MaskFilterType::New();
		maskFilter->SetInput(convertedToITK.outComposedImage->GetOutput());//input is a Image Pointer!!!
		maskFilter->SetMaskThreshold(m_MaskThreshold / m_DicomHelper->scaleSlope);//Get from UI or user-interaction
		maskFilter->Update();//output is a Image Pointer!!!
		//std::cout << "maskFilter: vectorLength = " << maskFilter->GetOutput()->GetVectorLength() << std::endl;

		//DTI claculation
		cout << "DTI is calculating.." << endl;
		typedef itk::GetDiffusionImageFilter<float, float> GetDiffusionImageType;
		GetDiffusionImageType::Pointer GetDiffusionImages = GetDiffusionImageType::New();
		GetDiffusionImages->SetInput(maskFilter->GetOutput());
		GetDiffusionImages->SetHMatrix(m_DicomHelper->finalH);
		GetDiffusionImages->SetSlice2PatMatrix(m_DicomHelper->slice2PatMatrix);
		GetDiffusionImages->SetBValueList(m_DicomHelper->BvalueList);
		GetDiffusionImages->Update();

		///////////////////////////////////////////
		//ITK to VTK
		typedef itk::ImageToVTKImageFilter <itk::Image<itk::RGBPixel< unsigned char >, 3>> itkToVtkConverter;
		itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
		convItkToVtk->SetInput(GetDiffusionImages->GetOutput2());
		convItkToVtk->Update();

		imageData->DeepCopy(convItkToVtk->GetOutput());
		cout << "color FA is finished calculating" << endl;
}
void DiffusionCore::EAdcCalculator(vtkSmartPointer <vtkImageData> imageData)
{
	/////////////////////////////////////
	//preprocessing: from VTK to ITK vector images
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
	//std::cout << "---------------------------- VOI is correct ? ---------------------" << std::endl;
	//this->DisplayDicomInfo(ExtractVOI->GetOutput());
	//std::cout << "---------------------------- End of VOI ---------------------" << std::endl;

	//PR: mask image doesn't change for slices other than the first one whlie sliding the maskthreshold slider
	//fix extent range here, default z extent is m_SourceImageCurrentSlice to m_SourceImageCurrentSlice;
	//Update image extent, otherwise maskFilter won't work for slices other than the first one (because seed starts at origin)
	//Update the origin as well
	vtkSmartPointer <vtkImageChangeInformation> changeInfo = vtkSmartPointer <vtkImageChangeInformation>::New();
	changeInfo->SetInputData(ExtractVOI->GetOutput());
	changeInfo->SetOutputOrigin(0, 0, 0);
	changeInfo->SetExtentTranslation(0, 0, -m_SourceImageCurrentSlice);
	changeInfo->Update();
	//std::cout << "---------------------------- translateExtent is correct ? ---------------------" << std::endl;
	//this->DisplayDicomInfo(changeInfo->GetOutput());
	//std::cout << "---------------------------- End of translateExtent ---------------------" << std::endl;



	if (m_DicomHelper->tensorComputationPossible)
	{
		//convert data from VTK to ITK: transfer to float and scaling to None scaled data.
		vtkSmartPointer <vtkImageData> imageSliceData = vtkSmartPointer <vtkImageData>::New();
		imageSliceData = changeInfo->GetOutput();
		VTKtoITK convertedToITK(m_DicomHelper->IsoImageLabel,
			m_DicomHelper->scaleSlope, m_DicomHelper->scaleIntercept, m_DicomHelper->numberOfComponents);
		convertedToITK.Scaling(imageSliceData);

		//Mask image filter
		typedef itk::MaskVectorImageFilter <float> MaskFilterType;
		MaskFilterType::Pointer maskFilter = MaskFilterType::New();
		maskFilter->SetInput(convertedToITK.outComposedImage->GetOutput());//input is a Image Pointer!!!
		maskFilter->SetMaskThreshold(m_MaskThreshold / m_DicomHelper->scaleSlope);//Get from UI or user-interaction
		maskFilter->Update();//output is a Image Pointer!!!
		//std::cout << "maskFilter: vectorLength = " << maskFilter->GetOutput()->GetVectorLength() << std::endl;

		//DTI claculation
		cout << "DTI is calculating.." << endl;
		typedef itk::GetDiffusionImageFilter<float, float> GetDiffusionImageType;
		GetDiffusionImageType::Pointer GetDiffusionImages = GetDiffusionImageType::New();
		GetDiffusionImages->SetInput(maskFilter->GetOutput());
		GetDiffusionImages->SetHMatrix(m_DicomHelper->finalH);
		GetDiffusionImages->SetSlice2PatMatrix(m_DicomHelper->slice2PatMatrix);
		GetDiffusionImages->SetBValueList(m_DicomHelper->BvalueList);
		GetDiffusionImages->Update();

		//Adc
		typedef itk::VectorIndexSelectionCastImageFilter <itk::VectorImage<float, 3>, FloatImageType> VectorImageToImageType;
		VectorImageToImageType::Pointer vectorImageToImageFilterAdc = VectorImageToImageType::New();
		vectorImageToImageFilterAdc->SetIndex(1);
		vectorImageToImageFilterAdc->SetInput(GetDiffusionImages->GetOutput1());
		vectorImageToImageFilterAdc->Update();

		//Rescale signal intensity to display
		typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
		RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
		rescaleFilter->SetInput(vectorImageToImageFilterAdc->GetOutput());
		rescaleFilter->SetOutputMaximum(4095.0);
		rescaleFilter->SetOutputMinimum(0.0);
		rescaleFilter->Update();
		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		//Data Clipping
		//typedef itk::DisplayOptimizer < FloatImageType, SrcImageType> DisplayOptimizerType;
		//DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
		//displayOptimizer->SetInput(rescaleFilter->GetOutput());
		//displayOptimizer->SetCoveragePercent(0.9);//Default is 0.99
		//displayOptimizer->Update();

		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		///////////////////////////////////////////
		//ITK to VTK
		typedef itk::ImageToVTKImageFilter <FloatImageType> itkToVtkConverter;
		itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
		convItkToVtk->SetInput(rescaleFilter->GetOutput());
		convItkToVtk->Update();

		imageData->DeepCopy(convItkToVtk->GetOutput());

	}
	else
	{
		for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i++)
		{
			//Handle each scalar component indivisually
			vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
			scalarComponent->SetInputData(changeInfo->GetOutput());
			scalarComponent->SetComponents(i);
			scalarComponent->Update();//Crutial, otherwise abort after running

			//VTK to ITK Image Data
			VtkToItkConverterType::Pointer vtkToItkImageFilter = VtkToItkConverterType::New();
			vtkToItkImageFilter->SetInput(scalarComponent->GetOutput());
			vtkToItkImageFilter->Update();

			//unsigned short image to float image
			CastFilterType::Pointer castFilter = CastFilterType::New();
			castFilter->SetInput(vtkToItkImageFilter->GetOutput());
			castFilter->Update();

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

		//get vector Image
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

		//std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
		// Adc Map filter
		typedef itk::AdcMapFilter <float, float> AdcMapFilterType;
		AdcMapFilterType::Pointer adcMap = AdcMapFilterType::New();
		adcMap->SetInput(maskFilter->GetOutput());
		adcMap->SetBValueList(this->m_DicomHelper->BvalueList);
		adcMap->Update();

		//std::cout << "------------- AdcMapFilter end runing ----------- " << std::endl;
		//std::cout << "adcVectorImage: vectorLength = " << adcMap->GetOutput()->GetVectorLength() << std::endl;

		//vector image to scalar image or imageContainer
		typedef itk::VectorIndexSelectionCastImageFilter <itk::VectorImage<float, 3>, FloatImageType> VectorImageToImageType;
		VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
		vectorImageToImageFilter->SetIndex(1);
		vectorImageToImageFilter->SetInput(adcMap->GetOutput());
		vectorImageToImageFilter->Update();

		//Rescale signal intensity to display
		typedef itk::RescaleIntensityImageFilter < FloatImageType, FloatImageType> RescaleIntensityImageType;
		RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
		rescaleFilter->SetInput(vectorImageToImageFilter->GetOutput());
		rescaleFilter->SetOutputMaximum(4095.0);
		rescaleFilter->SetOutputMinimum(0.0);
		rescaleFilter->Update();
		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		//Data Clipping
		typedef itk::DisplayOptimizer < FloatImageType, SrcImageType> DisplayOptimizerType;
		DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
		displayOptimizer->SetInput(rescaleFilter->GetOutput());
		displayOptimizer->SetCoveragePercent(0.9);//Default is 0.99
		displayOptimizer->Update();

		//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
		//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

		///////////////////////////////////////////
		//ITK to VTK
		typedef itk::ImageToVTKImageFilter <SrcImageType> itkToVtkConverter;
		itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
		convItkToVtk->SetInput(displayOptimizer->GetOutput());
		convItkToVtk->Update();

		imageData->DeepCopy(convItkToVtk->GetOutput());
		//Crutial not imageData = convItkToVtk->GetOutput(),
		//Because convItkToVtk->Update() pointer is recycled after this caller!!!

		//Visualize data
		//QVTKWidget *vtkComputedDwiWindow = new QVTKWidget;
		//this->m_Controls->displayLayout_2->addWidget(vtkComputedDwiWindow, 0, 1, 1, 1);
		//this->m_Controls->displayLayout_2->update();
		//QuantitativeImageViewer2D(convItkToVtk->GetOutput(), vtkComputedDwiWindow);
	}


}