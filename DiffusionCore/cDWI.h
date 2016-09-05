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

//DicomReader not DicomImageReader
#include "vtkDICOMDirectory.h"
#include "vtkDICOMItem.h"
#include "vtkDICOMMetaData.h"
#include "vtkDICOMDictionary.h"
#include "vtkDICOMReader.h"
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkIntArray.h>
#include "vnl_matrix.h"
#include "vnl_matrix_inverse.h"
#include "vnl_vector.h"
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
#include "itkAdcMapFilter.h"
#include "itkComputedDwiFilter.h"
#include "itkDwiIVIMFilter.h"


//#include "QuickView.h"
//#include "itkImageFileWriter.h"
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
	//int _Contrast;//ajust image window & level failed

public:
	void SetImageViewer(vtkImageViewer2* imageViewer) {
		_ImageViewer = imageViewer;
		_MinSlice = imageViewer->GetSliceMin();
		_MaxSlice = imageViewer->GetSliceMax();
		_Slice = _MinSlice;
		cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << std::endl;
	}

	void SetStatusMapper(vtkTextMapper* statusMapper) {
		_StatusMapper = statusMapper;
	}


protected:
	void MoveSliceForward() {
		if (_Slice < _MaxSlice) {
			_Slice += 1;
			cout << "MoveSliceForward::Slice = " << _Slice << std::endl;
			_ImageViewer->SetSlice(_Slice);
			std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
			_StatusMapper->SetInput(msg.c_str());
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
			_ImageViewer->Render();
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
		// forward event
		vtkInteractorStyleImage::OnKeyDown();
	}


	virtual void OnMouseWheelForward() {
		//std::cout << "Scrolled mouse wheel forward." << std::endl;
		MoveSliceForward();
		// don't forward events, otherwise the image will be zoomed 
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		// vtkInteractorStyleImage::OnMouseWheelForward();
	}


	virtual void OnMouseWheelBackward() {
		//std::cout << "Scrolled mouse wheel backward." << std::endl;
		if (_Slice > _MinSlice) {
			MoveSliceBackward();
		}
		// don't forward events, otherwise the image will be zoomed 
		// in case another interactorstyle is used (e.g. trackballstyle, ...)
		// vtkInteractorStyleImage::OnMouseWheelBackward();
	}

	/*virtual void OnLeftButtonDown(){
	_Contrast = 1;
	}

	virtual void OnLeftButtonUp(){
	_Contrast = 0;
	}*/

	virtual void OnMouseMove(){
		//if (_Contrast)
		/*{
		vtkSmartPointer < vtkRenderWindowInteractor> interactor = vtkSmartPointer < vtkRenderWindowInteractor>::New();
		interactor = _ImageViewer->GetRenderWindow()->GetInteractor();
		int lastPos[2];
		interactor->GetLastEventPosition(lastPos);
		int currPos[2];
		interactor->GetEventPosition(currPos);
		int deltaX = abs(lastPos[0] - currPos[0]);
		int deltaY = abs(lastPos[1] - currPos[1]);
		//reslice -> Update();
		double XSpacing = _ImageViewer->GetInput()->GetSpacing()[0];
		double YSpacing = _ImageViewer->GetInput()->GetSpacing()[1];
		cout << "Xsapcing: " << XSpacing*deltaX << endl;
		_ImageViewer->SetColorLevel(deltaX*XSpacing);
		_ImageViewer->SetColorWindow(deltaY*YSpacing);
		//}else*/
		{
			vtkInteractorStyle::OnMouseMove();
		}
	}

};

vtkStandardNewMacro(myVtkInteractorStyleImage);

class vtkSliderWindowCallback : public vtkCommand
{
public:
	static vtkSliderWindowCallback *New()
	{
		return new vtkSliderWindowCallback;
	}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkSliderWidget *sliderWidgetWindow =
			reinterpret_cast<vtkSliderWidget*>(caller);
		int valueWindow = static_cast<int>(static_cast<vtkSliderRepresentation *>(sliderWidgetWindow->GetRepresentation())->GetValue());
		this->_ImageViewer->SetColorWindow(valueWindow);
	}
	vtkSliderWindowCallback() :_ImageViewer(0) {}
	vtkImageViewer2* _ImageViewer;
};

class vtkSliderLevelCallback : public vtkCommand
{
public:
	static vtkSliderLevelCallback *New()
	{
		return new vtkSliderLevelCallback;
	}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkSliderWidget *sliderWidgetLevel =
			reinterpret_cast<vtkSliderWidget*>(caller);
		int valueLevel = static_cast<int>(static_cast<vtkSliderRepresentation *>(sliderWidgetLevel->GetRepresentation())->GetValue());
		this->_ImageViewer->SetColorLevel(valueLevel);
	}
	vtkSliderLevelCallback() :_ImageViewer(0) {}
	vtkImageViewer2* _ImageViewer;
};

#endif