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
//#include "vtkSmartPointer.h"
//#include "vtkImageTracerWidget.h"
//#include "vtkRendererCollection.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkPNGWriter.h"
#include "QVTKWidget.h"
#include "QVTKInteractor.h"

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
//#include <VTKtoITK.h> // not needed now
#include "vtkImageAppendComponents.h"
#include "vtkScalarBarActor.h"
#include <QMouseEvent>


//TRY SetNumberOfThreads(1) to solve the multi-thread ranmdom results
//Notify Wenxing

DiffusionCore::DiffusionCore(QWidget *parent)
	:QWidget(parent)
	//, dicomHelp(nullptr)
{	
	this->m_Controls = nullptr;
	this->sourceImage = vtkSmartPointer<vtkImageData>::New();//VTK image pointer
	this->m_MaskVectorImage = DiffusionCalculatorVectorImageType::New();

	m_SourceImageCurrentSlice = 0;
	m_QuantitativeImageCurrentSlice = 0;
	
	CreateQtPartControl(this);

	this->m_DicomHelper = nullptr;
	this->m_MaskThreshold = 3;
	this->m_ComputedBValue = 2000;
}

DiffusionCore::~DiffusionCore()
{
	//delete m_Controls;
	//this->m_Controls = NULL;

	this->sourceImage->Delete();
	this->sourceImage = NULL;
	this->m_MaskVectorImage->Delete();
	this->m_MaskVectorImage = ITK_NULLPTR;
}

void DiffusionCore::CreateQtPartControl(QWidget *parent)
{
	if (!m_Controls)
	{
		this->m_Controls = new Ui::DiffusionModule;
		this->m_Controls->setupUi(parent);
		m_Controls->bSlider->setMaximum(10000); //maximum supported b value;
		m_Controls->bSlider->setMinimum(500);
		m_Controls->bSlider->setDecimals(0);
		m_Controls->bSlider->setSingleStep(100);
		m_Controls->bSlider->setTickInterval(100);
		m_Controls->bSlider->setValue(2000);
		m_Controls->bSlider->setTracking(false);
		m_Controls->ThreshSlider->setMaximum(100); //maximum threshhold value;
		m_Controls->ThreshSlider->setDecimals(0);
		m_Controls->ThreshSlider->setSingleStep(1);
		m_Controls->ThreshSlider->setTickInterval(1);
		m_Controls->ThreshSlider->setValue(10);
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

		//Connect Toggles
		connect(m_Controls->adcToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcADC(bool)));
		connect(m_Controls->eadcToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcEADC(bool)));
		connect(m_Controls->cdwiToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcCDWI(bool)));
		connect(m_Controls->faToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcFA(bool)));
		connect(m_Controls->colorFAToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcColorFA(bool)));
		connect(m_Controls->ivimToggle, SIGNAL(toggled(bool)), this, SLOT(onCalcIVIM(bool)));

		//Connect Sliders
		connect(m_Controls->ThreshSlider, SIGNAL(valueChanged(double)), this, SLOT(onThreshSlide(double)));
		connect(m_Controls->bSlider, SIGNAL(valueChanged(double)), this, SLOT(onBSlide(double)));
		
		//Connect Pointer Buttons
		connect(m_Controls->pointer1, SIGNAL(toggled(bool)), this, SLOT(onRoiPointer(bool)));

		m_Controls->cDWI->hide();
		m_Controls->Thresh->hide();		
	}
}

void DiffusionCore::OnImageFilesLoaded(const QStringList& fileLists)
{	
	//enable all the buttons 
	this->m_Controls->ADCTool->setEnabled(true);
	this->m_Controls->DTITool->setEnabled(true);

	//Check if Original Images exist or not
	//QLayoutItem *existItem = this->m_Controls->displayLayout->itemAtPosition(0, 0);
	//QWidget * existWidget = existItem->widget();
	std::cout << this->m_Controls->displayLayout->rowCount() << "-" << this->m_Controls->displayLayout->columnCount() << std::endl;
	//if (this->m_Controls->displayLayout->rowCount() > 0 || this->m_Controls->displayLayout->columnCount() > 0) //if original image exists, delete all image windows. 
	if (layoutTable.size()>0)
	{
		for (int i = 0; i < layoutTable.size(); i++)
		{
			for (int j = 0; j < layoutTable[i].size(); j++)
			{
				ui_dumpWindow(i, j);
			}
		}	

		layoutTable.erase(layoutTable.begin(), layoutTable.end());				

		m_Controls->adcToggle->blockSignals(true);
		m_Controls->adcToggle->setChecked(false);
		m_Controls->adcToggle->blockSignals(false);

		m_Controls->eadcToggle->blockSignals(true);
		m_Controls->eadcToggle->setChecked(false);
		m_Controls->eadcToggle->blockSignals(false);

		m_Controls->cdwiToggle->blockSignals(true);
		m_Controls->cdwiToggle->setChecked(false);
		m_Controls->cdwiToggle->blockSignals(false);

		m_Controls->faToggle->blockSignals(true);
		m_Controls->faToggle->setChecked(false);
		m_Controls->faToggle->blockSignals(false);

		m_Controls->colorFAToggle->blockSignals(true);
		m_Controls->colorFAToggle->setChecked(false);
		m_Controls->colorFAToggle->blockSignals(false);

		m_Controls->ivimToggle->blockSignals(true);
		m_Controls->ivimToggle->setChecked(false);
		m_Controls->ivimToggle->blockSignals(false);

		this->m_Controls->adcToggle->setEnabled(true);
		this->m_Controls->eadcToggle->setEnabled(true);
		this->m_Controls->cdwiToggle->setEnabled(true);
		this->m_Controls->faToggle->setEnabled(true);
		this->m_Controls->colorFAToggle->setEnabled(true);
		this->m_Controls->ivimToggle->setEnabled(true);

		m_SourceImageCurrentSlice = 0;
		m_QuantitativeImageCurrentSlice = 0;
		this->m_DicomHelper = nullptr;
		this->m_MaskThreshold = 3;
		this->m_ComputedBValue = 2000;
	}

	//load data
	vtkStringArray* loadingFiles = vtkStringArray::New();
	loadingFiles->SetNumberOfValues(fileLists.size());
	for (int i = 0; i < fileLists.size(); i++)
	{
		loadingFiles->SetValue(i, fileLists.at(i).toStdString().c_str());
	}

	this->m_DicomHelper = new DicomHelper(loadingFiles);
	this->sourceImage = m_DicomHelper->DicomReader->GetOutput();
	
	this->UpdateMaskVectorImage();

	//tmp Test here to check if thresh slider is hided.
	this->m_Controls->Thresh->setVisible(false);
 
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
		}

		if (m_DicomHelper->numberOfBValue < 4)
		{
			this->m_Controls->ivimToggle->setDisabled(true);
		}
	}
	
	QVTKWidget *vtkWindow1 = new QVTKWidget;
	this->m_Controls->displayLayout->addWidget(vtkWindow1, 0, 0);
	//update the layout table
	std::vector<int> row1; row1.push_back(0);
	layoutTable.push_back(row1);

	std::vector<std::vector<int>>::iterator rowit;
	std::vector<int>::iterator colit;
	for (rowit = layoutTable.begin(); rowit != layoutTable.end(); rowit++)
	{
		std::cout << "[GRIDLAYOUT DEBUG INFO] BEFORE REORDER WINDOW>>>>>>>> |";
		for (colit = (*rowit).begin(); colit != (*rowit).end(); ++colit)
		{
			std::cout << *colit << " ";
		}
		std::cout << "|" << std::endl;
	}

	emit SignalDicomLoaded(true);

	std::cout << "srcimage viewer" << std::endl;
	DisplayDicomInfo(this->sourceImage);
	SourceImageViewer2D(this->sourceImage, vtkWindow1);
}

void DiffusionCore::UpdateMaskVectorImage()
{
	if (!this->sourceImage)
	{
		std::cout << "sourceImage NULL" << std::endl;
		return;
	}

	typedef itk::VectorContainer< SourceImagePixelType, DiffusionCalculatorImageType::Pointer > ImageContainerType;
	typedef itk::VTKImageToImageFilter <SourceImageType>	VtkToItkConverterType;
	typedef itk::CastImageFilter< SourceImageType, DiffusionCalculatorImageType >	CastFilterType;
	typedef itk::ShiftScaleImageFilter<DiffusionCalculatorImageType, DiffusionCalculatorImageType>	ShiftScaleType;


	ImageContainerType::Pointer imageContainer = ImageContainerType::New();
	if (m_DicomHelper->tensorComputationPossible && (m_DicomHelper->IsoImageLabel > -1))
	{		
		imageContainer->Reserve(this->m_DicomHelper->numberOfComponents -1);// only works for DTI b value 2, otherwise wrong
	}
	else
	{
		imageContainer->Reserve(this->m_DicomHelper->numberOfComponents);
	}
	//std::cout << "m_DicomHelper numof components" << this->m_DicomHelper->numberOfComponents << std::endl;

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
	//std::cout << " update vector image DTI numberOfComponents " << m_DicomHelper->numberOfComponents << std::endl;

	int DTIindex = 0;
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
		//DTI
		if (m_DicomHelper->tensorComputationPossible )
		{
			//std::cout << " update vector image DTI isoimagelabel " << m_DicomHelper->IsoImageLabel << std::endl;

			if (m_DicomHelper->IsoImageLabel != i)//Only handles numOfBvalue = 2, if > 2, it's wrong!!!
			{
				//imageContainer->Reserve(this->m_DicomHelper->numberOfComponents);
				imageContainer->InsertElement(DTIindex, dynamic_cast <DiffusionCalculatorImageType*> (shiftScale->GetOutput()));
				DTIindex++;
			}
		}
		else
		{
			//std::cout << " Update vector image DWI" << std::endl;
			imageContainer->InsertElement(i, dynamic_cast <DiffusionCalculatorImageType*> (shiftScale->GetOutput()));
		}		
	}

	//std::cout << " Update vector image DTI00  imageContainer size " << imageContainer->Size() << std::endl;
	//Get vector Image
	typedef itk::ComposeImageFilter<DiffusionCalculatorImageType>		ImageToVectorImageType;
	ImageToVectorImageType::Pointer imageToVectorImageFilter = ImageToVectorImageType::New();

	for (int i = 0; i < imageContainer->Size(); i++)
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
	typedef itk::MaskVectorImageFilter <DiffusionCalculatorPixelType> MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetInput(imageToVectorImageFilter->GetOutput());//input is a Image Pointer!!!
	maskFilter->SetMaskThreshold(m_MaskThreshold);//Get from UI or user-interaction
	maskFilter->Update();//output is a Image Pointer!!!

	//std::cout << "maskFilter: vectorLength = " << maskFilter->GetOutput()->GetVectorLength() << std::endl;

	//itk version of DeepCopy	
	this->m_MaskVectorImage->SetSpacing(maskFilter->GetOutput()->GetSpacing());
	this->m_MaskVectorImage->SetOrigin(maskFilter->GetOutput()->GetOrigin());
	this->m_MaskVectorImage->SetDirection(maskFilter->GetOutput()->GetDirection());
	this->m_MaskVectorImage->SetRegions(maskFilter->GetOutput()->GetLargestPossibleRegion());
	this->m_MaskVectorImage->SetVectorLength(maskFilter->GetOutput()->GetVectorLength());
	//std::cout << " m_MaskVectorImage length 00 0000000= " << this->m_MaskVectorImage->GetVectorLength() << std::endl;
	this->m_MaskVectorImage->Allocate();
	//std::cout << " m_MaskVectorImage length 00 = " << this->m_MaskVectorImage->GetVectorLength() << std::endl;


	typedef itk::ImageRegionConstIterator <DiffusionCalculatorVectorImageType> ConstMaskFilterIteratorType;
	typedef itk::ImageRegionIterator <DiffusionCalculatorVectorImageType> MaskVectorImageIteratorType;

	ConstMaskFilterIteratorType maskFilterIterator(maskFilter->GetOutput(), maskFilter->GetOutput()->GetLargestPossibleRegion());
	MaskVectorImageIteratorType maskVectorImageIterator(this->m_MaskVectorImage, this->m_MaskVectorImage->GetLargestPossibleRegion());

	maskFilterIterator.GoToBegin();
	maskVectorImageIterator.GoToBegin();
	while (!maskFilterIterator.IsAtEnd())
	{
		maskVectorImageIterator.Set(maskFilterIterator.Get());
		++maskFilterIterator;
		++maskVectorImageIterator;
	}
	//std::cout << " m_MaskVectorImage length 11 = " << this->m_MaskVectorImage->GetVectorLength() << std::endl;
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
		//this->m_Controls->Thresh->setVisible(_istoggled);

		//initialize first: clear up the widget 
		//Or should we just clear the renderer in current window?

		ui_RemoveWindow(ADC);
		ShareWindowEvent();

	}
	
}

void DiffusionCore::AdcCalculator(vtkSmartPointer <vtkImageData> imageData)
{	
	if (!this->m_MaskVectorImage) return;	
	//std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
	// Adc Map filter
	typedef itk::AdcMapFilter <DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> AdcMapFilterType;
	AdcMapFilterType::Pointer adcMap = AdcMapFilterType::New();
	adcMap->SetInput(this->m_MaskVectorImage);
	adcMap->SetBValueList(this->m_DicomHelper->BvalueList);
	adcMap->Update();
	//std::cout << "------------- AdcMapFilter end runing ----------- " << std::endl;
	//std::cout << "adcVectorImage: vectorLength = " << adcMap->GetOutput()->GetVectorLength() << std::endl;

	//vector image to scalar image or imageContainer
	typedef itk::VectorIndexSelectionCastImageFilter <DiffusionCalculatorVectorImageType, DiffusionCalculatorImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
	vectorImageToImageFilter->SetIndex(1);
	vectorImageToImageFilter->SetInput(adcMap->GetOutput());
	vectorImageToImageFilter->Update();


	//Rescale signal intensity to display
	typedef itk::RescaleIntensityImageFilter < DiffusionCalculatorImageType, DiffusionCalculatorImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
	rescaleFilter->SetInput(vectorImageToImageFilter->GetOutput());
	rescaleFilter->SetOutputMaximum(4095.0);
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->Update();

	m_ScalingParameter[2*ADC] = 1 / rescaleFilter->GetScale();
	m_ScalingParameter[2*ADC+1] = -1 * rescaleFilter->GetShift() / rescaleFilter->GetScale();
	std::cout << "rescaleFilter: inputMaximum = " << m_ScalingParameter[2 * ADC] << std::endl;
	std::cout << "rescaleFilter: inputMinimum = " << m_ScalingParameter[2 * ADC + 1] << std::endl;

	//Data Clipping
	typedef itk::DisplayOptimizer < DiffusionCalculatorImageType, SourceImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
	displayOptimizer->SetInput(rescaleFilter->GetOutput());
	displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer->Update();

	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	///////////////////////////////////////////
	//ITK to VTK
	typedef itk::ImageToVTKImageFilter <SourceImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
	convItkToVtk->SetInput(displayOptimizer->GetOutput());
	convItkToVtk->Update();

	imageData->DeepCopy(convItkToVtk->GetOutput());
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
		//this->m_Controls->Thresh->setVisible(_istoggled);
		
		//clear up existing widget
		ui_RemoveWindow(EADC);
		ShareWindowEvent();
	}
}

void DiffusionCore::EAdcCalculator(vtkSmartPointer <vtkImageData> imageData)
{
	if (!this->m_MaskVectorImage) return;
	//std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
	// Adc Map filter
	typedef itk::AdcMapFilter <DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> AdcMapFilterType;
	AdcMapFilterType::Pointer adcMap = AdcMapFilterType::New();
	adcMap->SetInput(this->m_MaskVectorImage);
	adcMap->SetBValueList(this->m_DicomHelper->BvalueList);
	adcMap->Update();
	//std::cout << "------------- AdcMapFilter end runing ----------- " << std::endl;
	//std::cout << "adcVectorImage: vectorLength = " << adcMap->GetOutput()->GetVectorLength() << std::endl;

	typedef itk::ComputedEadcFilter <DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> ComputedEadcFilterType;
	ComputedEadcFilterType::Pointer computedEadc = ComputedEadcFilterType::New();
	computedEadc->SetInput(adcMap->GetOutput());
	computedEadc->SetNumOfDiffDirections(this->m_DicomHelper->numberOfGradDirection);
	computedEadc->SetEadcBValue(this->m_DicomHelper->BvalueList.at(this->m_DicomHelper->BvalueList.size() - 1));//Get from UI input
	computedEadc->Update();

	//vector image to scalar image or imageContainer
	typedef itk::VectorIndexSelectionCastImageFilter <DiffusionCalculatorVectorImageType, DiffusionCalculatorImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
	vectorImageToImageFilter->SetIndex(0);
	vectorImageToImageFilter->SetInput(computedEadc->GetOutput());
	vectorImageToImageFilter->Update();


	//Rescale signal intensity to display
	typedef itk::RescaleIntensityImageFilter < DiffusionCalculatorImageType, DiffusionCalculatorImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
	rescaleFilter->SetInput(vectorImageToImageFilter->GetOutput());
	rescaleFilter->SetOutputMaximum(4095.0);
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->Update();
	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	//Data Clipping
	typedef itk::DisplayOptimizer < DiffusionCalculatorImageType, SourceImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
	displayOptimizer->SetInput(rescaleFilter->GetOutput());
	displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer->Update();

	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	///////////////////////////////////////////
	//ITK to VTK
	typedef itk::ImageToVTKImageFilter <SourceImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
	convItkToVtk->SetInput(displayOptimizer->GetOutput());
	convItkToVtk->Update();

	imageData->DeepCopy(convItkToVtk->GetOutput());
}

void DiffusionCore::onCalcCDWI(bool _istoggled)  
{


	if (_istoggled)
	{
		//Set Threshhold bar visible
		this->m_Controls->Thresh->setVisible(_istoggled);
		this->m_Controls->cDWI->setVisible(_istoggled);		

		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, CDWI);
		
		vtkSmartPointer <vtkImageData> computedDwi = vtkSmartPointer <vtkImageData>::New();
		this->CDWICalculator(computedDwi);
		QuantitativeImageViewer2D(computedDwi, vtkWindow,"cDWI");
	}
	else
	{
		//Hide Threshhold bar
		//this->m_Controls->Thresh->setVisible(_istoggled);
		this->m_Controls->cDWI->setVisible(_istoggled);		
		//clear up existing widget
		ui_RemoveWindow(CDWI);
		ShareWindowEvent();
	}
}

void DiffusionCore::CDWICalculator(vtkSmartPointer <vtkImageData> imageData)
{
	if (!this->m_MaskVectorImage) return;
	//std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
	// Adc Map filter
	typedef itk::AdcMapFilter <DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> AdcMapFilterType;
	AdcMapFilterType::Pointer adcMap = AdcMapFilterType::New();
	adcMap->SetInput(this->m_MaskVectorImage);
	adcMap->SetBValueList(this->m_DicomHelper->BvalueList);
	adcMap->Update();
	//std::cout << "------------- AdcMapFilter begin runing ----------- " << std::endl;
	//std::cout << "adcVectorImage: vectorLength = " << adcMap->GetOutput()->GetVectorLength() << std::endl;

	//cDwi filter
	typedef itk::ComputedDwiFilter <DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> ComputedDwiFilterType;
	ComputedDwiFilterType::Pointer computedDwi = ComputedDwiFilterType::New();
	computedDwi->SetInput(adcMap->GetOutput());
	computedDwi->SetNumOfDiffDirections(this->m_DicomHelper->numberOfGradDirection);
	computedDwi->SetComputedBValue(m_ComputedBValue);//Get from UI input
	computedDwi->Update();
	//std::cout << "cDWi vectorImage: vectorLength = " << computedDwi->GetOutput()->GetVectorLength() << std::endl;


	//vector image to scalar image or imageContainer
	typedef itk::VectorIndexSelectionCastImageFilter <DiffusionCalculatorVectorImageType, DiffusionCalculatorImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
	vectorImageToImageFilter->SetIndex(0);
	vectorImageToImageFilter->SetInput(computedDwi->GetOutput());
	vectorImageToImageFilter->Update();


	//Rescale signal intensity to display
	typedef itk::RescaleIntensityImageFilter < DiffusionCalculatorImageType, DiffusionCalculatorImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
	rescaleFilter->SetInput(vectorImageToImageFilter->GetOutput());
	rescaleFilter->SetOutputMaximum(4095.0);
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->Update();
	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	//Data Clipping
	typedef itk::DisplayOptimizer < DiffusionCalculatorImageType, SourceImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
	displayOptimizer->SetInput(rescaleFilter->GetOutput());
	displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer->Update();

	///////////////////////////////////////////
	//ITK to VTK for visualization
	typedef itk::ImageToVTKImageFilter < SourceImageType > itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
	convItkToVtk->SetInput(displayOptimizer->GetOutput());
	convItkToVtk->Update();

	imageData->DeepCopy(convItkToVtk->GetOutput());
}

void DiffusionCore::onThreshSlide(double maskThreshold) //It's a SLOT!!!
{
	m_MaskThreshold = maskThreshold/5;
	DiffusionCore::ShareWindowEvent();
}

void DiffusionCore::onBSlide(double computedBValue) //It's a SLOT!!!
{
	m_ComputedBValue = computedBValue;
	int cdwiRow(0), cdwiCol(0);
	ui_findWdw(CDWI, cdwiRow, cdwiCol);
	
	std::cout << "CDWI is at " << cdwiRow << ":" << cdwiCol << std::endl;

	//clear up existing widget
	QLayoutItem *existItem = this->m_Controls->displayLayout->itemAtPosition(cdwiRow, cdwiCol);
	QWidget * existWidget = existItem->widget();
	if (existWidget != NULL) 
	{
		QVTKWidget *cDwiWidget= static_cast <QVTKWidget*> (existWidget);
		vtkSmartPointer <vtkImageData> computedDwi = vtkSmartPointer <vtkImageData>::New();
		this->CDWICalculator(computedDwi);
		QuantitativeImageViewer2D(computedDwi, cDwiWidget, "cDWI");
	}
}

void DiffusionCore::onCalcFA(bool _istoggled)
{
	if (_istoggled)
	{
		//Set Threshhold bar visible
		this->m_Controls->Thresh->setVisible(_istoggled);

		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, FA);
		
		vtkSmartPointer <vtkImageData> calculatedFA = vtkSmartPointer <vtkImageData>::New();
		this->FaCalculator(calculatedFA);

		QuantitativeImageViewer2D(calculatedFA, vtkWindow, "FA");

	}
	else
	{
		//Hide Threshhold bar
		//this->m_Controls->Thresh->setVisible(_istoggled);

		ui_RemoveWindow(FA);
		ShareWindowEvent();
	}
}

void DiffusionCore::FaCalculator(vtkSmartPointer <vtkImageData> imageData)
{
	//this->UpdateMaskVectorImage();
	if (!this->m_MaskVectorImage) return;

	//DTI claculation
	cout << "DTI is calculating.." << endl;
	typedef itk::GetDiffusionImageFilter<DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> GetDiffusionImageType;
	GetDiffusionImageType::Pointer GetDiffusionImages = GetDiffusionImageType::New();
	GetDiffusionImages->SetInput(this->m_MaskVectorImage);
	GetDiffusionImages->SetHMatrix(m_DicomHelper->finalH);
	GetDiffusionImages->SetSlice2PatMatrix(m_DicomHelper->slice2PatMatrix);
	GetDiffusionImages->SetBValueList(m_DicomHelper->BvalueList);
	GetDiffusionImages->Update();

	//Adc
	typedef itk::VectorIndexSelectionCastImageFilter <DiffusionCalculatorVectorImageType, DiffusionCalculatorImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilterAdc = VectorImageToImageType::New();
	vectorImageToImageFilterAdc->SetIndex(2);
	vectorImageToImageFilterAdc->SetInput(GetDiffusionImages->GetOutput1());
	vectorImageToImageFilterAdc->Update();

	//Rescale signal intensity to display
	typedef itk::RescaleIntensityImageFilter < DiffusionCalculatorImageType, DiffusionCalculatorImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
	rescaleFilter->SetInput(vectorImageToImageFilterAdc->GetOutput());
	rescaleFilter->SetOutputMaximum(4095.0);
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->Update();
	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	//Data Clipping
	typedef itk::DisplayOptimizer < DiffusionCalculatorImageType, SourceImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
	displayOptimizer->SetInput(rescaleFilter->GetOutput());
	displayOptimizer->SetCoveragePercent(0.9);//Default is 0.99
	displayOptimizer->Update();

	//std::cout << "rescaleFilter: inputMaximum = " << rescaleFilter->GetInputMaximum() << std::endl;
	//std::cout << "rescaleFilter: inputMinimum = " << rescaleFilter->GetInputMinimum() << std::endl;

	///////////////////////////////////////////
	//ITK to VTK
	typedef itk::ImageToVTKImageFilter <SourceImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
	convItkToVtk->SetInput(displayOptimizer->GetOutput());
	convItkToVtk->Update();

	imageData->DeepCopy(convItkToVtk->GetOutput());
}

void DiffusionCore::onCalcColorFA(bool _istoggled)
{
	if (_istoggled)
	{
		//Set Threshhold bar visible
		this->m_Controls->Thresh->setVisible(_istoggled);

		QVTKWidget *vtkWindow = new QVTKWidget;

		int rowInd(0), colInd(0);
		ui_InsertWindow(rowInd, colInd, vtkWindow, CFA);

		vtkSmartPointer <vtkImageData> calculatedColorFA = vtkSmartPointer <vtkImageData>::New();
		this->ColorFACalculator(calculatedColorFA);
		QuantitativeImageViewer2D(calculatedColorFA, vtkWindow, "Color FA" );
	}
	else
	{
		//Hide Threshhold bar
		//this->m_Controls->Thresh->setVisible(_istoggled);

		ui_RemoveWindow(CFA);
		ShareWindowEvent();
	}
}

void DiffusionCore::ColorFACalculator(vtkSmartPointer <vtkImageData> imageData)
{
	//this->UpdateMaskVectorImage();
	if (!this->m_MaskVectorImage) return;

	//DTI claculation
	cout << "DTI is calculating.." << endl;
	typedef itk::GetDiffusionImageFilter<DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> GetDiffusionImageType;
	GetDiffusionImageType::Pointer GetDiffusionImages = GetDiffusionImageType::New();
	GetDiffusionImages->SetInput(this->m_MaskVectorImage);
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

void DiffusionCore::onCalcIVIM(bool _istoggled)
{
	if (_istoggled)
	{
		if (m_DicomHelper->BvalueList.size() < 4)
		{
			QMessageBox::StandardButton reply = QMessageBox::information(this, tr("QMessageBox::information()"), tr("Diffusion b valules less than 4, failed to calculate IVIM"));
			this->m_Controls->ivimToggle->blockSignals(true);
			this->m_Controls->ivimToggle->setChecked(false);
			this->m_Controls->ivimToggle->blockSignals(false);
			return;
		}

		if (m_DicomHelper->numberOfGradDirection > 1) // replace with extract the last direction data in future
		{
			QMessageBox::StandardButton reply = QMessageBox::information(this, tr("QMessageBox::information()"), tr("Diffusion Directions must be 1, failed to calculate IVIM"));
			this->m_Controls->ivimToggle->blockSignals(true);
			this->m_Controls->ivimToggle->setChecked(false);
			this->m_Controls->ivimToggle->blockSignals(false);
			return;
		}

		if (m_DicomHelper->BvalueList.at(m_DicomHelper->BvalueList.size() - 1) < 200) // replace with extract the last direction data in future
		{
			QMessageBox::StandardButton reply = QMessageBox::information(this, tr("QMessageBox::information()"), tr("Max b value must be larger than 200 s/mm2, failed to calculate IVIM"));
			this->m_Controls->ivimToggle->blockSignals(true);
			this->m_Controls->ivimToggle->setChecked(false);
			this->m_Controls->ivimToggle->blockSignals(false);
			return;
		}

		m_Controls->adcToggle->blockSignals(true);
		m_Controls->adcToggle->setChecked(false);
		m_Controls->adcToggle->blockSignals(false);

		m_Controls->eadcToggle->blockSignals(true);
		m_Controls->eadcToggle->setChecked(false);
		m_Controls->eadcToggle->blockSignals(false);

		m_Controls->cdwiToggle->blockSignals(true);
		m_Controls->cdwiToggle->setChecked(false);
		m_Controls->cdwiToggle->blockSignals(false);

		m_Controls->faToggle->blockSignals(true);
		m_Controls->faToggle->setChecked(false);
		m_Controls->faToggle->blockSignals(false);

		m_Controls->colorFAToggle->blockSignals(true);
		m_Controls->colorFAToggle->setChecked(false);
		m_Controls->colorFAToggle->blockSignals(false);

		this->m_Controls->adcToggle->setEnabled(false);
		this->m_Controls->eadcToggle->setEnabled(false);
		this->m_Controls->cdwiToggle->setEnabled(false);
		this->m_Controls->faToggle->setEnabled(false);
		this->m_Controls->colorFAToggle->setEnabled(false);
			
		vtkSmartPointer <vtkImageData> computedIVIM = vtkSmartPointer <vtkImageData>::New();
		this->IVIMCalculator(computedIVIM);
		if (!computedIVIM)
		{
			return;
		}
		
		int rowInd(0), colInd(0);
		QVTKWidget *vtkWindow = new QVTKWidget;
		ui_InsertWindow(rowInd, colInd, vtkWindow, IVIM_F); 
		IVIMImageViewer(computedIVIM, vtkWindow, 0); //IVIM_F = 0

		QVTKWidget *vtkWindow3 = new QVTKWidget;
		ui_InsertWindow(rowInd, colInd, vtkWindow3, IVIM_Dstar);
		IVIMImageViewer(computedIVIM, vtkWindow3, 1);  //IVIM_Dstar = 0

		QVTKWidget *vtkWindow2 = new QVTKWidget;
		ui_InsertWindow(rowInd, colInd, vtkWindow2, IVIM_D);
		IVIMImageViewer(computedIVIM, vtkWindow2, 2);  //IVIM_D = 0
	}
	else
	{ 
		ui_RemoveWindow(IVIM_Dstar);
		ui_RemoveWindow(IVIM_D);
		ui_RemoveWindow(IVIM_F);
		ShareWindowEvent();
	}
}

void DiffusionCore::IVIMCalculator(vtkSmartPointer <vtkImageData> imageData)
{
	if (!this->m_MaskVectorImage) return;

	//IVIM calculator
	std::cout << "------------- IVIM begin runing ----------- " << std::endl;
	typedef itk::DwiIVIMFilter2 <DiffusionCalculatorPixelType, DiffusionCalculatorPixelType> IVIMFilterType;
	IVIMFilterType::Pointer dwiIVIM = IVIMFilterType::New();
	dwiIVIM->SetInput(this->m_MaskVectorImage);
	dwiIVIM->SetBValueList(this->m_DicomHelper->BvalueList);
	dwiIVIM->SetNumOfIterations(50);
	dwiIVIM->Update();
	std::cout << "------------- IVIM end runing ----------- " << std::endl;
	std::cout << "IVIM: vectorLength = " << dwiIVIM->GetOutput()->GetVectorLength() << std::endl;

	//Component 0---------------------------------------
	//vector image to scalar image or imageContainer
	typedef itk::VectorIndexSelectionCastImageFilter <DiffusionCalculatorVectorImageType, DiffusionCalculatorImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilter = VectorImageToImageType::New();
	vectorImageToImageFilter->SetIndex(0);
	vectorImageToImageFilter->SetInput(dwiIVIM->GetOutput());
	vectorImageToImageFilter->Update();

	//Data Clipping
	typedef itk::DisplayOptimizer < DiffusionCalculatorImageType, DiffusionCalculatorImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer = DisplayOptimizerType::New();
	displayOptimizer->SetInput(vectorImageToImageFilter->GetOutput());
	displayOptimizer->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer->Update();

	//Rescale signal intensity to display
	//It doesn't take vector image type, so handle each component seperately
	typedef itk::RescaleIntensityImageFilter < DiffusionCalculatorImageType, DiffusionCalculatorImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter = RescaleIntensityImageType::New();
	rescaleFilter->SetInput(displayOptimizer->GetOutput());
	rescaleFilter->SetOutputMaximum(255.0);
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->Update();

	///////////////////////////////////////////
	//ITK to VTK for visualization
	typedef itk::ImageToVTKImageFilter < DiffusionCalculatorImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk = itkToVtkConverter::New();
	convItkToVtk->SetInput(rescaleFilter->GetOutput());// changed to rescale filter
	convItkToVtk->Update();
	//}

	vtkSmartPointer <vtkImageAppendComponents> appendComponents = vtkSmartPointer <vtkImageAppendComponents>::New();
	appendComponents->SetInputData(convItkToVtk->GetOutput());

	//Component 1---------------------------------------
	//vector image to scalar image or imageContainer
	//typedef itk::VectorIndexSelectionCastImageFilter <DiffusionCalculatorVectorImageType, DiffusionCalculatorImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilter1 = VectorImageToImageType::New();
	vectorImageToImageFilter1->SetIndex(1);
	vectorImageToImageFilter1->SetInput(dwiIVIM->GetOutput());
	vectorImageToImageFilter1->Update();

	//Data Clipping
	//typedef itk::DisplayOptimizer < DiffusionCalculatorImageType, DiffusionCalculatorImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer1 = DisplayOptimizerType::New();
	displayOptimizer1->SetInput(vectorImageToImageFilter1->GetOutput());
	displayOptimizer1->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer1->Update();

	//Rescale signal intensity to display
	//It doesn't take vector image type, so handle each component seperately
	//typedef itk::RescaleIntensityImageFilter < DiffusionCalculatorImageType, DiffusionCalculatorImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter1 = RescaleIntensityImageType::New();
	rescaleFilter1->SetInput(displayOptimizer1->GetOutput());
	rescaleFilter1->SetOutputMaximum(255.0);
	rescaleFilter1->SetOutputMinimum(0.0);
	rescaleFilter1->Update();

	///////////////////////////////////////////
	//ITK to VTK for visualization
	//typedef itk::ImageToVTKImageFilter < DiffusionCalculatorImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk1 = itkToVtkConverter::New();
	convItkToVtk1->SetInput(rescaleFilter1->GetOutput());// changed to rescale filter
	convItkToVtk1->Update();
	
	appendComponents->AddInputData(convItkToVtk1->GetOutput());

	//Component 2---------------------------------------
	//vector image to scalar image or imageContainer
	//typedef itk::VectorIndexSelectionCastImageFilter <DiffusionCalculatorVectorImageType, DiffusionCalculatorImageType> VectorImageToImageType;
	VectorImageToImageType::Pointer vectorImageToImageFilter2 = VectorImageToImageType::New();
	vectorImageToImageFilter2->SetIndex(2);
	vectorImageToImageFilter2->SetInput(dwiIVIM->GetOutput());
	vectorImageToImageFilter2->Update();

	//Data Clipping
	//typedef itk::DisplayOptimizer < DiffusionCalculatorImageType, DiffusionCalculatorImageType> DisplayOptimizerType;
	DisplayOptimizerType::Pointer displayOptimizer2 = DisplayOptimizerType::New();
	displayOptimizer2->SetInput(vectorImageToImageFilter2->GetOutput());
	displayOptimizer2->SetCoveragePercent(0.98);//Default is 0.99
	displayOptimizer2->Update();

	//Rescale signal intensity to display
	//It doesn't take vector image type, so handle each component seperately
	//typedef itk::RescaleIntensityImageFilter < DiffusionCalculatorImageType, DiffusionCalculatorImageType> RescaleIntensityImageType;
	RescaleIntensityImageType::Pointer rescaleFilter2 = RescaleIntensityImageType::New();
	rescaleFilter2->SetInput(displayOptimizer2->GetOutput());
	rescaleFilter2->SetOutputMaximum(255.0);
	rescaleFilter2->SetOutputMinimum(0.0);
	rescaleFilter2->Update();

	///////////////////////////////////////////
	//ITK to VTK for visualization
	//typedef itk::ImageToVTKImageFilter < DiffusionCalculatorImageType> itkToVtkConverter;
	itkToVtkConverter::Pointer convItkToVtk2 = itkToVtkConverter::New();
	convItkToVtk2->SetInput(rescaleFilter2->GetOutput());// changed to rescale filter
	convItkToVtk2->Update();

	appendComponents->AddInputData(convItkToVtk2->GetOutput());
	appendComponents->Update();
	std::cout << "IVIM calculator imageData components before deepcopy = " << std::endl;
	imageData->DeepCopy(appendComponents->GetOutput());
	std::cout << "IVIM calculator imageData components = " << imageData->GetNumberOfScalarComponents() << std::endl;
}

void DiffusionCore::onCursorPickValue(vtkObject* obj)
{

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

	this->m_Controls->displayLayout->addWidget(vtkWindow, rowInd, colInd);


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

	//this->m_Controls->displayLayout->update();
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

	//for (rowit = layoutTable.begin(),rowCounter=0; rowit != layoutTable.end(); rowit++, rowCounter++)
	//{
	//	for (colit = (*rowit).begin(),colCounter=0; colit != (*rowit).end();)
	//	{
	//		//std::cout << "Running at " << rowCounter << ":" << colCounter << std::endl;
	//		if (*colit == imageLabel)
	//		{
	//			colit = (*rowit).erase(colit);
	//			del_rowInd = rowCounter;
	//			del_colInd = colCounter;
	//		}
	//		else
	//		{
	//			++colit; colCounter++;
	//		}
	//	}
	//}	
	
	ui_findWdw(imageLabel, del_rowInd, del_colInd);	
	layoutTable[del_rowInd].erase(remove(layoutTable[del_rowInd].begin(), layoutTable[del_rowInd].end(), imageLabel), layoutTable[del_rowInd].end());
	std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>deleting the widget type "<< imageLabel <<" at "<<del_rowInd<<":"<<del_colInd << std:: endl;
	std::vector<int>::iterator rmIt = layoutTable[del_rowInd].begin() + del_colInd;


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

		if (del_rowInd == cur_RowNum - 1)//if it is the last row. 
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

	//this->m_Controls->displayLayout->update();	

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
	std::cout << "[GRIDLAYOUT DEBUG INFO]>>>>>>>>deleting window " << " at " << row << ":" << col << std::endl;
	QLayoutItem *existItem = this->m_Controls->displayLayout->itemAtPosition(row, col);
	QWidget * existWidget = existItem->widget();
	std::cout << "There is a widget ? " << (existWidget != NULL) <<std::endl;
	if (existWidget != NULL) {		
		this->m_Controls->displayLayout->removeWidget(existWidget);
		existWidget->setParent(NULL);//if you want to delete the widget, do: widget->setParent(NULL); delete widget; 
		delete existWidget;
	}
	//this->m_Controls->displayLayout->update();
	//std::cout << m_Controls->displayLayout->rowCount() << " Row " << m_Controls->displayLayout->columnCount() << " Col " << std::endl;
}

void DiffusionCore::ui_findWdw(imageType imageLabel, int& row, int& col)
{
	//std::vector<std::vector<int>>::iterator rowit;
	std::vector<int>::iterator colit;
	int rowCounter(0), colCounter(0);
	std::cout << "Find image " << imageLabel;
	//STEP1, Find window rendering assigned image and delete it, update the layoutTable. 
	for (rowCounter = 0; rowCounter < layoutTable.size(); rowCounter++)
	{
		for (colit = layoutTable[rowCounter].begin(), colCounter = 0; colit < layoutTable[rowCounter].end(); colit++, colCounter++)
		{

			if (*colit == imageLabel)
			{
				row = rowCounter;
				col = colCounter;
			}
		}
		std::cout << " row = " << row << "col = " << col << std::endl;
	}
}

void DiffusionCore::ShareWindowEvent()
{

	if (layoutTable[0].size() > 1)
	{
		this->UpdateMaskVectorImage();
		if (!this->m_MaskVectorImage) return;
	}

	std::cout << "[GRIDLAYOUT DEBUG INFO] Before Rendering girdlayout is ";
	std::cout << layoutTable.size() << " X " << layoutTable[0].size() << std::endl;

	vtkSmartPointer <vtkImageData> ivimImage = vtkSmartPointer <vtkImageData>::New();
	int row(0), col(0);
	ui_findWdw(IVIM_D, row, col);
	if (row != 0 || col != 0)
	{
		this->IVIMCalculator(ivimImage);
	}

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
					vtkSmartPointer <vtkImageData> thisImage = vtkSmartPointer <vtkImageData>::New();
					QVTKWidget *thisWindow = static_cast <QVTKWidget*> (existWidget);
					switch (layoutTable[rowCounter][colCounter])
					{
					case ORIGINAL:
						break;
					case ADC:
						this->AdcCalculator(thisImage);
						QuantitativeImageViewer2D(thisImage, thisWindow, "ADC");
						break;
					case CDWI:
						this->CDWICalculator(thisImage);
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
					case IVIM_F:						
						IVIMImageViewer(ivimImage, thisWindow, 0);
						break;
					case IVIM_Dstar:
						IVIMImageViewer(ivimImage, thisWindow, 1);
						break;
					case IVIM_D:
						IVIMImageViewer(ivimImage, thisWindow, 2);
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

void DiffusionCore::onRoiPointer(bool _istoggled)
{
	if (_istoggled)
	{
		//loop over all current QVTKWidgets
		for (int rowCounter = 0; rowCounter < layoutTable.size(); rowCounter++)
		{
			for (int colCounter = 0; colCounter < layoutTable[rowCounter].size(); colCounter++)
			{
				if (colCounter != 0 || rowCounter != 0)
				{
					QLayoutItem *existItem = this->m_Controls->displayLayout->itemAtPosition(rowCounter, colCounter);
					QWidget * existWidget = existItem->widget();
					QVTKWidget *thisWindow = static_cast <QVTKWidget*> (existWidget);

					thisWindow->setAutomaticImageCacheEnabled(true);

					vtkSmartPointer<vtkImageTracerWidget> tracer =
						vtkSmartPointer<vtkImageTracerWidget>::New();
					tracer->GetLineProperty()->SetRepresentationToSurface();
					tracer->SetCaptureRadius(30);
					tracer->GetLineProperty()->SetLineWidth(1.5);
					tracer->GetLineProperty()->SetColor(1, 0, 0);
					tracer->ProjectToPlaneOn();
					tracer->SetInteractor(thisWindow->GetRenderWindow()->GetInteractor());
					tracer->SetViewProp((vtkProp*)thisWindow->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetViewProps()->GetItemAsObject(0));
					tracer->AutoCloseOn();

					Connections = vtkEventQtSlotConnect::New();

					// get right mouse pressed with high priority

					//Connections->Connect(thisWindow->GetRenderWindow()->GetInteractor(), vtkCommand::EndInteractionEvent, 
					//	this, SLOT(onDisplayROIstat(thisWindow->cachedImage, tracer, layoutTable[rowCounter][colCounter])));

					Connections->Connect(tracer, vtkCommand::EndInteractionEvent,
						this, SLOT(onDisplayROIstat(vtkObject*, unsigned long, void*, void*, vtkCommand*)), thisWindow);

					tracer->On();
				}
			}
		}
	}
	else
	{
		ShareWindowEvent();
	}

}


void DiffusionCore::onDisplayROIstat(vtkObject* tracer, unsigned long, void* client_data, void*, vtkCommand* )
{
	vtkImageTracerWidget* _tracer = static_cast<vtkImageTracerWidget*>(tracer);

	QVTKWidget* _window = static_cast<QVTKWidget*>(client_data);
	//vtkWindowToImageFilter* imageGrabber = vtkWindowToImageFilter::New();
	//imageGrabber->SetInput(window->GetRenderWindow());
	//imageGrabber->Update();
	vtkImageActor* _actor = static_cast<vtkImageActor*>(_window->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetViewProps()->GetItemAsObject(0));

	vtkImageData* o_image = _actor->GetMapper()->GetInput();

	std::cout << "image dimensions are " << o_image->GetDimensions()[0] <<"-" <<o_image->GetDimensions()[1] << std::endl;
	std::cout << "image dimensions have components " << o_image->GetNumberOfScalarComponents() << std::endl;

	vtkSmartPointer<vtkImageExtractComponents> extractRedFilter =
		vtkSmartPointer<vtkImageExtractComponents>::New();
	extractRedFilter->SetInputData(o_image);
	extractRedFilter->SetComponents(0);
	extractRedFilter->Update();

	vtkImageData* _image = extractRedFilter->GetOutput();

	vtkSmartPointer<vtkPNGWriter> writer =
		vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName("E:\\screenshot2.png");
	writer->SetInputData(_image);
	writer->Write();
	
	if (!_tracer) { return; }

	if (!_tracer->IsClosed())
	{
		QString warnInfo(tr("ROI contour is not closed. Please Draw Again"));
		this->m_Controls->statisBrowser->setText(warnInfo);
		return;
	}

	vtkSmartPointer<vtkPolyData> path =
		vtkSmartPointer<vtkPolyData>::New();
	_tracer->GetPath(path);
	std::cout << "There are " << path->GetNumberOfPoints() << " points in the path." << std::endl;
	//vtkSmartPointer<vtkImageShiftScale> scale =
	//	vtkSmartPointer<vtkImageShiftScale>::New();
	//scale->SetInputData(imageData);
	//scale->SetScale(scalingPara[0]);
	//scale->SetShift(scalingPara[1]);
	//scale->Update();

	vtkSmartPointer<vtkPolyDataToImageStencil> polyDataToImageStencil =
		vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	polyDataToImageStencil->SetTolerance(0);
	polyDataToImageStencil->SetInputData(path); // if version < 5 setinputConnection
	polyDataToImageStencil->SetOutputOrigin(_image->GetOrigin());
	polyDataToImageStencil->SetOutputSpacing(_image->GetSpacing());
	polyDataToImageStencil->SetOutputWholeExtent(_image->GetExtent());
	polyDataToImageStencil->Update();

	vtkSmartPointer<vtkImageAccumulate> imageAccumulate =
		vtkSmartPointer<vtkImageAccumulate>::New();
	imageAccumulate->SetStencilData(polyDataToImageStencil->GetOutput());
	imageAccumulate->SetInputData(_image);
	imageAccumulate->Update();

	//calculate related parameters:
	float scalingValue, shiftValue;
	//switch (_imageLabel)
	//{
	//case ADC:
	//	scalingValue = scalingPara[0];
	//	shiftValue = scalingPara[1];
	//	break;
	//case FA:
	//	scalingValue = scalingPara[6];
	//	shiftValue = scalingPara[7];
	//	break;
	//}
	scalingValue = m_ScalingParameter[2*ADC];
	shiftValue = m_ScalingParameter[2 * ADC + 1];
	
	std::cout << "onROIdrawing: scalingValue = " << m_ScalingParameter[2 * ADC] << std::endl;
	std::cout << "onROIdrawing: shiftValue = " << m_ScalingParameter[2 * ADC + 1] << std::endl;

	float Area = imageAccumulate->GetVoxelCount()*_image->GetSpacing()[0] * _image->GetSpacing()[1];
	float Mean = *imageAccumulate->GetMean()*scalingValue + shiftValue;	
	float Max = *imageAccumulate->GetMax()*scalingValue + shiftValue;
	float Min = *imageAccumulate->GetMin()*scalingValue + shiftValue;
	float std = *imageAccumulate->GetStandardDeviation()*scalingValue + shiftValue;

	//cout << "Area:" << Area << endl;
	//cout << "Mean:" << Mean << endl;
	//cout << "Max:" << Max << endl;
	//cout << "Min:" << Min << endl;
	//cout << "std:" << std << endl;
	
	QString RoiInfo(tr("Area : %1").arg(Area));
	RoiInfo.append("\n");
	RoiInfo.append(tr("Mean : %1").arg(Mean));
	RoiInfo.append("\n");
	RoiInfo.append(tr("Max : %1").arg(Max));
	RoiInfo.append("\n");
	RoiInfo.append(tr("Min : %1").arg(Min));
	RoiInfo.append("\n");
	RoiInfo.append(tr("std : %1").arg(std));
	RoiInfo.append("\n");
	this->m_Controls->statisBrowser->setText(RoiInfo);
}


void DiffusionCore::DisplayDicomInfo(vtkSmartPointer <vtkImageData> imageData)
{

	//cout << "fileNames: " << reader->GetFileNames()<< endl;
	//qDebug() << "number of COMPONENTS: " << imageData->GetNumberOfScalarComponents() << endl;

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
	//qDebug() << "image dims: " << dims[0] << "x" << dims[1] << "x" << dims[2] << endl;
	imageData->GetOrigin(origins);
	//qDebug() << "image origins: " << origins[0] << " " << origins[1] << " " << origins[2] << endl;
	imageData->GetSpacing(spacing);
	//qDebug() << "image spacing: " << spacing[0] << "x" << spacing[1] << "x" << spacing[2] << endl;
	imageData->GetExtent(extent);
	//qDebug() << "extent: " << extent[0] << "x" << extent[1] << "x" << extent[2] << "x" << extent[3] << "x" << extent[4] << "x" << extent[5] << endl;
	imageData->GetScalarRange(range);//1. cannot be type of float here, it's a bug of vtk?  2. error while calculating quantitative output images: pipeline should be updated before calling this method!!!
	//qDebug() << "range: " << range[0] << "x" << range[1] << endl;
	//std::cout << " imageData: GetScalarType() = " << imageData->GetScalarType() << std::endl;
	//std::cout << " scaleSlope = " << this->m_DicomHelper->scaleSlope << std::endl;
	//std::cout << "scaleIntercept = " << this->m_DicomHelper->scaleIntercept << std::endl;

	//std::cout << "diffusion related parameters---" << std::endl;
	//for (int i = 0; i < this->m_DicomHelper->numberOfComponents; i = i + this->m_DicomHelper->numberOfGradDirection)
	//{
	//	int j = 0;
	//	std::cout << "bValueList " << j << ": " << this->m_DicomHelper->BvalueList.at(i / this->m_DicomHelper->numberOfGradDirection) << std::endl;
	//	j++;
	//}
	//std::cout << "numberOfComponents = " << this->m_DicomHelper->numberOfComponents << std::endl;
	//std::cout << "numberOfComponents = " << this->m_DicomHelper->numberOfComponents << std::endl;
	//std::cout << "numberOfGradDirection = " << this->m_DicomHelper->numberOfGradDirection << std::endl;
	//std::cout << "numberOfBValue = " << this->m_DicomHelper->numberOfBValue << std::endl;

	vtkMedicalImageProperties* properties = m_DicomHelper->DicomReader->GetMedicalImageProperties();
	QString imageInfo(tr("Patient Name : "));
	imageInfo.append(QLatin1String(properties->GetPatientName()));
	imageInfo.append("\n");
	imageInfo.append(tr("Study Name: "));
	imageInfo.append(QLatin1String(properties->GetStudyDescription()));
	imageInfo.append("\n");
	imageInfo.append(tr("Scan Name: "));
	imageInfo.append(QLatin1String(properties->GetSeriesDescription()));
	imageInfo.append("\n");
	imageInfo.append(tr("Scan Date: "));
	imageInfo.append(QLatin1String(properties->GetAcquisitionDate()));
	imageInfo.append("\n");

	imageInfo.append(tr("Image Dimension: [ %1 : %2 : %3 ]\n").arg(dims[0]).arg(dims[1]).arg(dims[2]));
	imageInfo.append(tr("Number Of b Value : %1 \n").arg(m_DicomHelper->numberOfBValue));
	imageInfo.append(tr("Number Of Gradient Direction : %1 \n").arg(m_DicomHelper->numberOfGradDirection));
	m_Controls->infoBrowser->setText(imageInfo);
}

void DiffusionCore::SourceImageViewer2D(vtkSmartPointer <vtkImageData> imageData, QVTKWidget *qvtkWidget)
{
	//double *imageDataRange = new double[2];
	//imageDataRange = imageData->GetScalarRange();//Replace with to be displayed

	vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
	imageViewer->SetInputData(imageData);
	imageViewer->SetSliceOrientationToXY();
	//imageViewer->SetColorWindow(imageDataRange[1] - imageDataRange[0]);
	//imageViewer->SetColorLevel(0.5* (imageDataRange[1] + imageDataRange[0]));

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

	//For Source imageviewer: use renderwindowInteractor. 
	//so that imageviewer pointer would not be released. Otherwise move slice forward/backward will report error, because imageviewer has been released.
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

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
	//std::cout << "RenderWindow SIZE = " << *(qvtkWidget->GetRenderWindow()->GetSize()) << std::endl;
	//std::cout << "ImageViewer SIZE = " << *(imageViewer->GetSize()) << std::endl;
	//std::cout << "QVTKWIDEGT SIZE = " << qvtkWidget->width() << "-" << qvtkWidget->height() << std::endl;
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetSize(800, 800);
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetPosition(qvtkWidget->x(), qvtkWidget->y());
	qvtkWidget->GetRenderWindow()->SetInteractor(renderWindowInteractor);//crutial to let qvtkWidget share the same interactor with imageViewer

	imageViewer->GetRenderer()->ResetCamera(); //Reset camera and then render is better
	vtkSmartPointer<vtkCamera> camera = imageViewer->GetRenderer()->GetActiveCamera();
	this->SetImageFillWindow(camera, imageData, qvtkWidget->width(), qvtkWidget->height());
	//imageViewer->GetRenderer()->SetActiveCamera(camera);
	//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
	qvtkWidget->GetRenderWindow()->Render();
	renderWindowInteractor->Initialize();
	qvtkWidget->show();
	renderWindowInteractor->Start();
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
	//std::cout << "QVTKWIDEGT SIZE = " << * (qvtkWidget->GetRenderWindow()->GetSize()) << std::endl;
	//std::cout << "ImageViewer SIZE = " << *(imageViewer->GetSize()) << std::endl;
	//std::cout << "QVTKWIDEGT SIZE = " << qvtkWidget->width()<< "-" << qvtkWidget->height() << std::endl;
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetSize(800, 800);
	//qvtkWidget->GetRenderWindow()->vtkRenderWindow::SetPosition(qvtkWidget->x(), qvtkWidget->y());
	qvtkWidget->GetRenderWindow()->SetInteractor(renderWindowInteractor);//crutial to let qvtkWidget share the same interactor with imageViewer
	//std::cout << "QVTKWIDEGT SIZE after interactor= " << qvtkWidget->width() << "-" << qvtkWidget->height() << std::endl;
	//imageViewer->SetupInteractor(renderWindowInteractor);

	imageViewer->GetRenderer()->ResetCamera(); //Reset camera and then render is better
	vtkSmartPointer<vtkCamera> camera = imageViewer->GetRenderer()->GetActiveCamera();
	this->SetImageFillWindow(camera, imageData, qvtkWidget->width(), qvtkWidget->height());
	//imageViewer->GetRenderer()->SetActiveCamera(camera);
	//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
	qvtkWidget->GetRenderWindow()->Render();
	renderWindowInteractor->Initialize();
	qvtkWidget->show();//qvtkWidget->show() changes window size!!!!
	//std::cout << "QVTKWIDEGT SIZE after show= " << qvtkWidget->width() << "-" << qvtkWidget->height() << std::endl;
}

void DiffusionCore::IVIMImageViewer(vtkSmartPointer <vtkImageData> imageData, QVTKWidget *qvtkWidget, int imageIdx)
{
	vtkSmartPointer <vtkLookupTable> lookupTable = vtkSmartPointer <vtkLookupTable>::New();
	lookupTable->SetNumberOfTableValues(256);//try below color range
	lookupTable->SetTableRange(0.1, 255.1);//try below color range
	lookupTable->SetHueRange(0.66667,0.0);//rainbow color map: from blue to red
	lookupTable->UseBelowRangeColorOn();
	lookupTable->SetBelowRangeColor(0.0,0.0,0.0,1.0);
	lookupTable->Build();

	//vtkSmartPointer<vtkColorTransferFunction> lookupTable = vtkSmartPointer<vtkColorTransferFunction>::New();
	//lookupTable->SetColorSpaceToRGB();
	//lookupTable->AddRGBPoint(0,0.0,0.0,1.0);
	//lookupTable->AddRGBPoint(255, 1.0, 0.0, 0.0);
	//lookupTable->SetScaleToLinear();
	vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
	scalarComponent->SetInputData(imageData);
	scalarComponent->SetComponents(imageIdx);
	scalarComponent->Update();

	vtkSmartPointer <vtkImageMapToColors> scalarValueToColors = vtkSmartPointer <vtkImageMapToColors>::New();
	scalarValueToColors->SetLookupTable(lookupTable);
	scalarValueToColors->PassAlphaToOutputOn();
	scalarValueToColors->SetInputData(scalarComponent->GetOutput());

	vtkSmartPointer <vtkImageActor> imageActor = vtkSmartPointer <vtkImageActor>::New();
	imageActor->GetMapper()->SetInputConnection(scalarValueToColors->GetOutputPort());
	imageActor->GetProperty()->SetInterpolationTypeToNearest();
	
	//imageActor->GetProperty()->SetAmbient(1);
	//imageActor->GetProperty()->SetDiffuse(0);

	vtkSmartPointer <vtkScalarBarActor> scalarBar = vtkSmartPointer <vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(lookupTable);
	//scalarBar->SetTitle("Title");
	//scalarBar->SetNumberOfLabels(4);//Default is 5
	scalarBar->SetDrawTickLabels(0);//Disable labels

	//Adjust scalarBar positon according to image Actor
	//double imageActorPos[3];
	//imageActor->GetPosition(imageActorPos);
	//double scalarBarPos[1];
	//scalarBar->GetPosition();
	//imageActor->GetMaxXBound();
	//std::cout << "imageActor pos = " << imageActorPos[0] << " " << imageActorPos[1] << " " << imageActorPos[2] << std::endl;
	//std::cout << "imageActor xBound = " << imageActor->GetMaxXBound() << " " << imageActor->GetMinXBound() << std::endl;
	//std::cout << "imageActor yBound = " << imageActor->GetMaxYBound() << " " << imageActor->GetMinYBound() << std::endl;

	//std::cout << "scalarBar default pos = " << scalarBar->GetPosition()[0] << " " << scalarBar->GetPosition()[1] << std::endl;
	//std::cout << "scalarBar default height width = " << scalarBar->GetHeight() << " " << scalarBar->GetWidth() << std::endl;
	//std::cout << "QVTKWIDEGT SIZE BEFORE show= " << qvtkWidget->width() << "-" << qvtkWidget->height() << std::endl;

	//std::cout << "imageActor xBound = " << imageActor->GetMaxXBound() << " " << imageActor->GetMinXBound() << std::endl;
	//imageActor->get//GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
	//scalarBar->GetPositionCoordinate()->SetCoordinateSystemToWorld();
	//scalarBar->GetPositionCoordinate()->SetValue(0.9 * imageActor->GetMaxXBound(), 0.8 * imageActor->GetMaxXBound());


	vtkSmartPointer <vtkRenderer> renderer = vtkSmartPointer <vtkRenderer>::New();
	renderer->AddActor2D(scalarBar);
	renderer->AddActor(imageActor);
	renderer->ResetCamera();
	

	//renderer->SetErase(1);
	//renderer->GradientBackgroundOn();
	//renderer->SetBackground2();
	renderer->SetBackground(0.0, 0.0, 0.0);

	vtkSmartPointer <vtkRenderWindow> renderWindow = vtkSmartPointer <vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	//renderWindow->SetBorders(1);
	//std::cout << "renderer size = " << renderer->GetSize()[0] << " " << renderer->GetSize()[1] << std::endl;
	vtkSmartPointer <QVTKInteractor> interactor = vtkSmartPointer <QVTKInteractor>::New();

	vtkSmartPointer <vtkInteractorStyleImage> interactorStyle = vtkSmartPointer <vtkInteractorStyleImage>::New();

	interactor->SetInteractorStyle(interactorStyle);

	qvtkWidget->SetRenderWindow(renderWindow);
	qvtkWidget->GetRenderWindow()->SetInteractor(interactor);


	renderer->ResetCamera();//Reset camera and then render is better
	vtkSmartPointer<vtkCamera> camera = renderer->GetActiveCamera();
	this->SetImageFillWindow(camera, scalarComponent->GetOutput(), qvtkWidget->width(), qvtkWidget->height());
	//imageViewer->GetRenderer()->SetActiveCamera(camera);
	//imageViewer->GetRenderer()->SetBackground(0.2, 0.3, 0.4);
	qvtkWidget->GetRenderWindow()->Render();
	interactor->Initialize();
	qvtkWidget->show();
	//std::cout << "QVTKWIDEGT SIZE AFTER show= " << qvtkWidget->width() << "-" << qvtkWidget->height() << std::endl;

}

void DiffusionCore::SetImageFillWindow(vtkSmartPointer <vtkCamera> & camera, vtkSmartPointer <vtkImageData> imageData, double width, double height){

	if (!(camera || imageData)) return;

	//int dims[3];
	double origins[3];
	double spacing[3];
	int extent[6];
//	double range[2];
	double center[3];
	//imageData->GetDimensions(dims);
	imageData->GetOrigin(origins);
	imageData->GetSpacing(spacing);
	imageData->GetExtent(extent);



	camera->ParallelProjectionOn();

	double xFov = (-extent[0] + extent[1] + 1)*spacing[0];
	double yFov = (-extent[2] + extent[3] + 1)*spacing[1];

	double screenRatio = height / width;
	double imageRatio = yFov / xFov;

	double parallelScale = imageRatio > screenRatio ? yFov : screenRatio/imageRatio*yFov;

	double focalDepth = camera->GetDistance();

	center[0] = origins[0] + spacing[0] * .5*(extent[0] + extent[1]);
	center[1] = origins[1] + spacing[1] * .5*(extent[2] + extent[3]);
	center[2] = origins[2] + spacing[2] * .5*(extent[4] + extent[5]);

	camera->SetParallelScale(0.5*parallelScale);
	//std::cout << "parallel scale minHalfFov " << minHalfFov << std::endl;
	//std::cout << "focal depth " << focalDepth << std::endl;
	//std::cout << "view angle " << camera->GetViewAngle() << std::endl;
	//std::cout << "Quantitative imageviewer, input width, input height =  " << width << " " << height << std::endl;
	//std::cout << "center " << center[0] << " " << center[1] << " " << center[2] << std::endl;
	camera->SetFocalPoint(center[0], center[1], center[2]);
	camera->SetPosition(center[0], center[1], focalDepth);
	//camera->SetViewAngle(80);
}