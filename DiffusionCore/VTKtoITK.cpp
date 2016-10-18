#include <VTKtoITK.h>
#include <vtkImageExtractComponents.h>


VTKtoITK::VTKtoITK(int &isolabel, float &scaleSlope, float &scaleInterp, int &componentsNumber)
{
	slope = scaleSlope;
	interp = scaleInterp;
	numberOfComponents = componentsNumber;
	isoImageLabel = isolabel;
	if (isoImageLabel > -1)
		vectorContainedNumber = numberOfComponents - 1;
	else
		vectorContainedNumber = numberOfComponents;

	imageContainer = ImageContainerType::New();
	imageContainer->Reserve(vectorContainedNumber);
 
	//cout << "VTK to ITK initialnizeing...." << endl;
};

void VTKtoITK::Scaling(vtkSmartPointer <vtkImageData> &output)
{
	int index = 0;
	vtkSmartPointer <vtkImageData> srcVTKData = vtkSmartPointer <vtkImageData>::New();
	srcVTKData = output;
	for (int i = 0; i < numberOfComponents; i++)
	{
		//Handle each scalar component indivisually
		vtkSmartPointer <vtkImageExtractComponents> scalarComponent = vtkSmartPointer <vtkImageExtractComponents>::New();
		scalarComponent->SetInputData(output);
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
		shiftScale->SetShift(-interp);
		shiftScale->SetScale(1.0 / slope);
		shiftScale->Update();

		//Save vector image & image container
		if (isoImageLabel != i)
		{		
			imageContainer->InsertElement(index, dynamic_cast <FloatImageType*> (shiftScale->GetOutput()));
			//outComposedImage->SetInput(index, dynamic_cast <FloatImageType*> (shiftScale->GetOutput()));
			++index;
		}
	}


	//cout << "image container calculated." << endl;
	//cout << "vector contained nubmer" << vectorContainedNumber << endl;
	for (int i = 0; i < vectorContainedNumber; i++)
	{
		outComposedImage->SetInput(i, imageContainer->GetElement(i));
	}
	outComposedImage->Update();
	//cout << "composed image calculated." << endl;
};


