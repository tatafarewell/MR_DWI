#include <itkImageToImageFilter.h>
#include <itkVTKImageToImageFilter.h>
#include <itkImage.h>
#include <itkVectorContainer.h>
#include <itkComposeImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkCastImageFilter.h>

#include <vtkAlgorithmOutput.h>
const int dimension = 3;
typedef itk::Image < float, dimension> FloatImageType;

class VTKtoITK
{

public:
	VTKtoITK::VTKtoITK(int &isolabel, float &scaleSlope, float &scaleInterp, int &componentsNumber);
	void VTKtoITK::Scaling(vtkSmartPointer <vtkImageData> &output);


	typedef itk::VectorContainer< unsigned int, FloatImageType::Pointer> ImageContainerType;
	ImageContainerType::Pointer imageContainer;
	typedef itk::ComposeImageFilter <FloatImageType> outComposedImageType;
	outComposedImageType::Pointer outComposedImage = outComposedImageType::New();

	typedef itk::Image < unsigned short, dimension> SrcImageType;
	typedef itk::VTKImageToImageFilter <SrcImageType>	VtkToItkConverterType;
	typedef itk::CastImageFilter< SrcImageType, FloatImageType >	CastFilterType;
	typedef itk::ShiftScaleImageFilter<FloatImageType, FloatImageType>	ShiftScaleType;

	
private:

	float slope;
	float interp;
	int numberOfComponents;
	int isoImageLabel;

	int vectorContainedNumber;
};

