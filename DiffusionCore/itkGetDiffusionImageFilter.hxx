#ifndef _itkGetDiffusionImageFilter_hxx_
#define _itkGetDiffusionImageFilter_hxx_

#include "itkGetDiffusionImageFilter.h"
#include "vnl\algo\vnl_matrix_inverse.h"

#include <itkAdcMapFilter.h>
#include <itkTensor.h>


namespace itk
{
	template< typename TInputPixelType, typename TOutputPixelType >
	GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >
		::GetDiffusionImageFilter()
	{
		this->SetNumberOfRequiredOutputs(3);
		this->SetNthOutput(0, this->MakeOutput(0));
		this->SetNthOutput(1, this->MakeOutput(1));
		this->SetNthOutput(2, this->MakeOutput(2));

		m_InputImage = ITK_NULLPTR;
		m_OutputImage = ITK_NULLPTR;
		m_OutputTensor = ITK_NULLPTR;
		m_OutputImageColorFa = ITK_NULLPTR;
		
	}

	template< typename TInputPixelType, typename TOutputPixelType >
	void
		GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >
		::BeforeThreadedGenerateData()
	{
		m_InputImage = this->GetInput();

		m_OutputImage = this->GetOutput1();
		m_OutputImage->SetVectorLength(3);
		m_OutputImage->Allocate();

		m_OutputImageColorFa = this->GetOutput2();

		m_OutputTensor = this->GetOutput3();
	}

	template< typename TInputPixelType, typename TOutputPixelType >
	void
		GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >
		::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId)
	{
		//calculate ADC map for all the gradient directions, note the ouput of this filter!
		typedef itk::AdcMapFilter<float, float> ADCFilterType;
		ADCFilterType::Pointer ADCMapFilter = ADCFilterType::New();
		ADCMapFilter->SetInput(m_InputImage);
		ADCMapFilter->SetBValueList(m_BValueList);
		ADCMapFilter->Update();

		InputImageType::Pointer ADCMap = ADCMapFilter->GetOutput();
		ConstInputIteratorType sourceImageIt(m_InputImage, outputRegionForThread);
		ConstInputIteratorType inputIt(ADCMap, outputRegionForThread);

		//OutputIteratorType outputItAdc(m_OutputImageAdc, outputRegionForThread);	
		OutputIteratorType outputIt(m_OutputImage, outputRegionForThread);
		OutputIteratorTypeColorFa outputItColorFa(m_OutputImageColorFa, outputRegionForThread);
		OutputIteratorTypeTensor outputTensorIt(m_OutputTensor, outputRegionForThread);

		colorFaPixelType colorFaImage;
		itk::Tensor<float> Tensors;

		vnl_matrix<float> combinationMatrix = vnl_matrix_inverse<float>(m_finalHMatrix.transpose()*m_finalHMatrix)
			*m_finalHMatrix.transpose();

		sourceImageIt.GoToBegin();
		inputIt.GoToBegin();
		//outputItAdc.GoToBegin();
		outputIt.GoToBegin();
		outputItColorFa.GoToBegin();

		float Adc, EAdc, Fa;
		colorFaPixelType colorFa;

		while (!inputIt.IsAtEnd())
		{
			float referenceValue = sourceImageIt.Get()[0];
			if (sourceImageIt.Get()[0] < exp(-10))
			{
				Adc = 0.0;
				Fa = 0.0;
				EAdc = 0.0;
				colorFa.SetRed(0.0);
				colorFa.SetBlue(0.0);
				colorFa.SetGreen(0.0);
				Tensors.Set(sourceImageIt.Get()[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
     		}
			else
			{
				int numberOfDirection = m_finalHMatrix.rows();
				vnl_vector<float> adcVector(numberOfDirection);
				vnl_vector<float> tensor(6);
				for (int i = 0; i < numberOfDirection; i++)
				{
					adcVector[i] = inputIt.Get()[2 * i + 1];
				}
				tensor = combinationMatrix*adcVector;
				
				
				Tensors.Set(sourceImageIt.Get()[0], tensor[0], tensor[1], tensor[2], tensor[3], tensor[4], tensor[5]);

				Adc = Tensors.ComputeAdcIso();
				Fa = Tensors.ComputeFa();
				EAdc = exp(-Adc*m_BValueList[1]);
					//Tensors.ComputeEAdcIso(m_BValueList[1]);
				colorFa = Tensors.ComputeColorFa(m_slice2PatMatrix);
			}

			VariableLengthVector<float> diffusionPara;
			diffusionPara.SetSize(3);
			diffusionPara[0] = Adc;
			diffusionPara[1] = EAdc;
			diffusionPara[2] = Fa;
			outputIt.Set(diffusionPara);
			++outputIt;

			outputItColorFa.Set(colorFa);
			++outputItColorFa;

			outputTensorIt.Set(Tensors);
			++outputTensorIt;

			++sourceImageIt;
			++inputIt;
		}
	}

	template< typename TInputPixelType, typename TOutputPixelType >
	DataObject::Pointer GetDiffusionImageFilter<TInputPixelType, TOutputPixelType>::MakeOutput(unsigned int idx)
	{
		DataObject::Pointer output;
		switch (idx)
		{
		case 0:
			output = (OutputImageType::New()).GetPointer();
			break;
		case 1:
			output = (colorFaImageType::New()).GetPointer();
			break;
		case 2:
			output = (tensorImageType::New()).GetPointer();
			break;
		default:
			std::cerr << "No output " << idx << std::endl;
			output = NULL;
			break;
		}
		return output.GetPointer();
	}

	template< typename TInputPixelType, typename TOutputPixelType >
	typename GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >::OutputImageType*
		GetDiffusionImageFilter <TInputPixelType, TOutputPixelType>::GetOutput1()
	{
		return dynamic_cast< OutputImageType * >(
			this->ProcessObject::GetOutput(0));
	}

	template< typename TInputPixelType, typename TOutputPixelType >
	typename GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >::colorFaImageType*
		GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >
	      ::GetOutput2()
	{
		return dynamic_cast< colorFaImageType * >(
			this->ProcessObject::GetOutput(1));
	}

	template< typename TInputPixelType, typename TOutputPixelType >
	typename GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >::tensorImageType*
		GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >
		::GetOutput3()
	{
		return dynamic_cast< tensorImageType * >(
			this->ProcessObject::GetOutput(2));
	}

	template< typename TInputPixelType, typename TOutputPixelType >
	void
		GetDiffusionImageFilter< TInputPixelType, TOutputPixelType >
		::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		//os << indent << "NumOfComponents: " << m_NumOfComponents << std::endl;
	}


} // end namespace itk

#endif