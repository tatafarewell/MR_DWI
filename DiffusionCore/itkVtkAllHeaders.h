#ifndef vtkcDWI_h
#define vtkcDWI_h

////////////////////////////////
//VTK headers: redundant~~
// some standard vtk headers
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
// headers needed for this example
#include <vtkImageViewer2.h>
#include <vtkDICOMImageReader.h>
#include <vtkInteractorStyleImage.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkTextMapper.h>

#include <vtkImageData.h>
#include <vtkImageMapper.h>
#include <vtkImageChangeInformation.h>
#include <vtkExtractVOI.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkImageActor.h>
#include <vtkOutlineFilter.h>
#include <vtkImageSliceMapper.h>
#include <vtkMatrix4x4.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkWarpScalar.h>
#include <vtkImagePlaneWidget.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkJPEGReader.h>
#include <vtkCommand.h>
#include <vtkSliderWidget.h>
#include <vtkWidgetEvent.h>
#include <vtkCallbackCommand.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProperty.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageShiftScale.h>
#include <vtkImageChangeInformation.h>
//jiangli test reset window
#include "vtkInformation.h"
#include "vtkAlgorithm.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkImageMapper3D.h"
#include "vtkImageMapToWindowLevelColors.h"
//jiangli test reset window end

//DicomReader not DicomImageReader
#include "vtkDICOMDirectory.h"
#include "vtkDICOMItem.h"
#include "vtkDICOMMetaData.h"
#include "vtkDICOMDictionary.h"
#include "vtkDICOMReader.h"
#include "vtkMedicalImageProperties.h"
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkIntArray.h>
#include <vtkImageProperty.h>
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_vector.h"

////////////////////////////////////////
//ITK headers
#include "itkImageToImageFilter.h"
#include "itkArray.h"
#include "itkSmartPointer.h"
#include "itkImage.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkVectorImage.h"
#include "itkVectorContainer.h"
#include "itkVariableLengthVector.h"
#include "itkComposeImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include <itkAdcMapFilter.h>
#include <itkComputedDwiFilter.h>
#include <itkDwiIVIMFilter.h>
#include <itkMaskVectorImageFilter.h>
#include <itkDisplayOptimizer.h>
#include <itkComputedEadcFilter.h>

#include "itkMatrix.h"
#include "itkArray.h"
//#include "itkNumericTraits.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"
//#include "itkImageDuplicator.h"
#include "itkShiftScaleImageFilter.h"
//ITK segmentation
#include "itkGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

////////////////////////////////////////
// Glue headers
// needed to easily convert from int to std::string
#include <sstream>
#include <string> //difference between string and string.h?
#include <stdlib.h>
#include <iostream>
#include <vector>


class StatusMessage {
public:
	static std::string Format(int slice, int maxSlice) {
		std::stringstream tmp;
		tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
		return tmp.str();
	}
};


// Define own interaction style
class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
	static myVtkInteractorStyleImage* New();
	vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);

protected:
	vtkImageViewer2* _ImageViewer;
	vtkTextMapper* _StatusMapper;
	int _Slice;
	int _MinSlice;
	int _MaxSlice;

	int *_CurrentSlice;//Change from outside

	int _Component;
	int _MaxComponent;

	vtkImageData* _OriginalInputImageData;

public:
	void SetImageViewer(vtkImageViewer2* imageViewer) {
		_ImageViewer = imageViewer;
		_MinSlice = imageViewer->GetSliceMin();
		_MaxSlice = imageViewer->GetSliceMax();
		_Slice = _MinSlice;

		_OriginalInputImageData = imageViewer->GetInput();
		
		SetDefaultWindowLevel(_Slice, _Component);

		_MaxComponent = imageViewer->GetInput()->GetNumberOfScalarComponents() - 1;
		_Component = 0;

		//cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << std::endl;
	}

	void SetStatusMapper(vtkTextMapper* statusMapper) {
		_StatusMapper = statusMapper;
	}

	void GetCurrentSliceNumber(int & CurrentSlice)
	{
		_CurrentSlice = &CurrentSlice;
	}

	void SetDefaultWindowLevel(int currentSlice, int currentComponent)
	{
		int dims[3];
		_OriginalInputImageData->GetDimensions(dims);
		if (dims[2] > 1)
		{
			//Extract slice
			vtkSmartPointer <vtkExtractVOI> ExtractVOI = vtkSmartPointer <vtkExtractVOI>::New();
			ExtractVOI->SetInputData(_OriginalInputImageData);
			ExtractVOI->SetVOI(0, dims[0] - 1, 0, dims[1] - 1, currentSlice, currentSlice);
			ExtractVOI->Update();

			//Extract component
			vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
			scalarComponent->SetInputData(ExtractVOI->GetOutput());
			scalarComponent->SetComponents(currentComponent);
			scalarComponent->Update();

			vtkSmartPointer <vtkImageChangeInformation> changeInfo = vtkSmartPointer <vtkImageChangeInformation>::New();
			changeInfo->SetInputData(scalarComponent->GetOutput());
			changeInfo->SetOutputOrigin(0, 0, 0);
			changeInfo->SetExtentTranslation(0, 0, -currentSlice);
			changeInfo->Update();
			std::cout << "before window leve range = " << std::endl;
			double *range = static_cast<double *>(changeInfo->GetOutput()->GetScalarRange());
			std::cout << "window leve range = " << range[0] << " " << range[1] << std::endl;
			_ImageViewer->SetColorWindow(range[1] - range[0]);
			_ImageViewer->SetColorLevel(0.5* (range[1] + range[0]));
		}
		else if (dims[2] == 1)
		{
			// move quantitative image viewer window level to here later
		}
	}

protected:
	void MoveSliceForward() {
		if (_Slice < _MaxSlice) {
			_Slice += 1;
			cout << "MoveSliceForward::Slice = " << _Slice << std::endl;

			_ImageViewer->SetSlice(_Slice);
			std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
			_StatusMapper->SetInput(msg.c_str());
			*_CurrentSlice = _Slice;

			SetDefaultWindowLevel(_Slice, _Component);

			_ImageViewer->Render();
		}
	}

	void MoveSliceBackward() {
		if (_Slice > _MinSlice) {
			_Slice -= 1;
			cout << "MoveSliceBackward::Slice = " << _Slice << std::endl;
			_ImageViewer->SetSlice(_Slice);
			std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
			_StatusMapper->SetInput(msg.c_str());
			*_CurrentSlice = _Slice;

			SetDefaultWindowLevel(_Slice, _Component);

			_ImageViewer->Render();			
		}
	}

	void MoveSliceComponentForward()
	{
		if (_Component < _MaxComponent)
		{
			_Component ++;
			cout << "MoveSliceComponentForward: Component = " << _Component << std::endl;

			vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
			scalarComponent->SetInputData(_OriginalInputImageData);
			scalarComponent->SetComponents(_Component);
			scalarComponent->Update();
			_ImageViewer->SetInputData(scalarComponent->GetOutput());
			_ImageViewer->SetSlice(_Slice);

			SetDefaultWindowLevel(_Slice, _Component);

			_ImageViewer->Render();

			//Needed to set input back???
			//_ImageViewer->SetInputData(_OriginalInputImageData);
		}
	}

	void MoveSliceComponentBackward()
	{
		if (_Component > 0)
		{
			_Component--;
			cout << "MoveSliceComponentBackward: Component = " << _Component << std::endl;
			vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
			scalarComponent->SetInputData(_OriginalInputImageData);
			scalarComponent->SetComponents(_Component);
			scalarComponent->Update();
			_ImageViewer->SetInputData(scalarComponent->GetOutput());
			_ImageViewer->SetSlice(_Slice);

			SetDefaultWindowLevel(_Slice, _Component);

			_ImageViewer->Render();

			//Needed to set input back???
			//_ImageViewer->SetInputData(_OriginalInputImageData);
		}
	}



	virtual void OnKeyDown() {
		std::string key = this->GetInteractor()->GetKeySym();
		if (key.compare("Up") == 0) {
			//cout << "Up arrow key was pressed." << endl;
			MoveSliceForward();
		}
		else if (key.compare("Down") == 0) {
			//cout << "Down arrow key was pressed." << endl;
			MoveSliceBackward();
		}
		else if (key.compare("Left") == 0){
			MoveSliceComponentBackward();
		}
		else if (key.compare("Right") == 0){
			MoveSliceComponentForward();
		}

		// forward event
		vtkInteractorStyleImage::OnKeyDown();
	}


	virtual void OnMouseWheelForward() {
		MoveSliceForward();
		// don't forward events, otherwise the image will be zoomed 
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		//vtkInteractorStyleImage::OnMouseWheelForward();
	}


	virtual void OnMouseWheelBackward() {
		//std::cout << "Scrolled mouse wheel backward." << std::endl;
		if (_Slice > _MinSlice) {
			MoveSliceBackward();
		}
		// don't forward events, otherwise the image will be zoomed 
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		//vtkInteractorStyleImage::OnMouseWheelBackward();
	}

	//Get keybord input, only support Reset window level for now
	virtual void OnChar()
	{
		switch (this->GetInteractor()->GetKeyCode())
		{
			case 'r':
			case 'R':
			{	
				//PR fix: press reset window on quantitative window results crash but ont on source image window
				//Because Source image use vtkRenderWindowInteractor type, Quantitative image uses QVTKInteractor type!!!
				//Plus, _ImageViewer->GetInput()->GetScalarRange() doesn't work for source image, why?
				if (_OriginalInputImageData)
				{
					//std::cout << "reset window 11" << std::endl;
					//std::cout << "image data quantitative viewer " << _ImageViewer->GetInput()->GetNumberOfScalarComponents() << std::endl;
					//std::cout << "image data quantitative viewer range 1 " << _OriginalInputImageData->GetScalarRange()[1] << std::endl;
					double *range = static_cast<double *>(_OriginalInputImageData->GetScalarRange());
					range[1] = (range[1] - range[0] > 255 ? 255 + range[0] : range[1]);
					// *range is from 0 to rescaleFilter->GetOutputMaximum().
					//Howerver, windows will automaticcaly change the display range to [0,255], thus if we reset color window to rescaleFilter->GetOutputMaximum(), it won't be optimal
					vtkImageProperty *property = this->CurrentImageProperty;
					property->SetColorWindow(range[1] - range[0]);					
					property->SetColorLevel(0.5 * (range[1] + range[0]));
					this->GetInteractor()->GetRenderWindow()->Render();
				}
				break;
			}
			default:
				//this->Superclass::OnChar();
				break;
		}
	}

	virtual void WindowLevel()
	{
		vtkRenderWindowInteractor *rwi = this->Interactor;

		this->WindowLevelCurrentPosition[0] = rwi->GetEventPosition()[0];
		this->WindowLevelCurrentPosition[1] = rwi->GetEventPosition()[1];

		if (this->HandleObservers &&
			this->HasObserver(vtkCommand::WindowLevelEvent))
		{
			this->InvokeEvent(vtkCommand::WindowLevelEvent, this);
		}
		else if (this->CurrentImageProperty)
		{
			int *size = this->CurrentRenderer->GetSize();

			double window = this->CurrentImageProperty->GetColorWindow();
			//double window = this->WindowLevelInitial[0];//WindwoLevelInitial doesn't work properly, replace it
			double level = this->CurrentImageProperty->GetColorLevel();
			// Compute normalized delta
			double dx = (this->WindowLevelCurrentPosition[0] -
				this->WindowLevelStartPosition[0]) * 0.00005;//jiangli modify 
			double dy = (this->WindowLevelStartPosition[1] -
				this->WindowLevelCurrentPosition[1]) * 0.00005; // / size[1];//jiangli modify

			// Scale by current values
			if (fabs(window) > 0.01)
			{
				dx = dx * window;
			}
			else
			{
				dx = dx * (window < 0 ? -0.01 : 0.01);
			}
			if (fabs(level) > 0.01)
			{
				dy = dy * level;
			}
			else
			{
				dy = dy * (level < 0 ? -0.01 : 0.01);
			}

			// Abs so that direction does not flip

			if (window < 0.0)
			{
				dx = -1 * dx;
			}
			if (level < 0.0)
			{
				dy = -1 * dy;
			}

			// Compute new window level

			double newWindow = dx + window;
			double newLevel = level - dy;

			if (newWindow < 0.01)
			{
				newWindow = 0.01;
			}
			newLevel = newLevel < 0.01 ? 0.01 : newLevel;

			this->CurrentImageProperty->SetColorWindow(newWindow);
			this->CurrentImageProperty->SetColorLevel(newLevel);
			this->Interactor->Render();
		}
	}
};

vtkStandardNewMacro(myVtkInteractorStyleImage);

class vtkTestCallbackCommand : public vtkCommand
{
public:

	static vtkTestCallbackCommand *New()
	{
		return new vtkTestCallbackCommand;
	}

	virtual void Execute(vtkObject *caller, unsigned long eventId, void*)
	{
		vtkSmartPointer <vtkRenderWindowInteractor> iren = static_cast <vtkRenderWindowInteractor*> (caller);

		std::string key = iren->GetKeySym();
		if (key.compare("Left") == 0)
		{
			std::cout << "TestCallbackFunc Left " << std::endl;
		}
		else if (key.compare("Right") == 0)
		{
			std::cout << "TestCallbackFunc Right " << std::endl;
		}
		else if (key.compare("Down") == 0)
		{
			std::cout << "TestCallbackFunc Down" << std::endl;
		}

		//int valueWindow = static_cast<int>(static_cast<vtkSliderRepresentation *>(sliderWidgetWindow->GetRepresentation())->GetValue());
		//this->_ImageViewer->SetColorWindow(valueWindow);
	}

	vtkTestCallbackCommand()
	{

	}
	//vtkImageViewer2* _ImageViewer;
};

#endif
