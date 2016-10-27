#include <DicomHelper.h>

#include <vtkDICOMMetaData.h>
#include <vtkDICOMValue.h>
#include <vtkImageData.h>
#include <vtkStringArray.h>

void DicomHelper::SetDiffusionSeriesNumber(int seriesNumber)
{
	DiffusionSeryNumber = seriesNumber;
};

void DicomHelper::DicomInfo()
{	
	DicomReader->Update();
	numberOfComponents = DicomReader->GetOutput()->GetNumberOfScalarComponents();
	DicomReader->GetOutput()->GetDimensions(imageDimensions);

	DicomReader->UpdateInformation();
	vtkSmartPointer<vtkDICOMMetaData> meta = vtkSmartPointer<vtkDICOMMetaData>::New();
	meta = DicomReader->GetMetaData();

	if (meta->HasAttribute( DICOM_SCALE_SLOP ))
	{
		scaleSlope = meta->GetAttributeValue(DICOM_SCALE_SLOP).AsFloat();
	}

	if (meta->HasAttribute( DICOM_SCALE_INTERPCEPT ))
	{
		scaleIntercept = meta->GetAttributeValue(DICOM_SCALE_INTERPCEPT).AsFloat();
	}

	if (meta->HasAttribute(DICOM_NO_DIFF_BVAL))
	{
	    numberOfBValue= meta->GetAttributeValue(DICOM_NO_DIFF_BVAL).AsInt();
	}

	if (meta->HasAttribute(DICOM_NO_DIFF_GRAD_ORIENT))
	{
		numberOfGradDirection = meta->GetAttributeValue(DICOM_NO_DIFF_GRAD_ORIENT).AsInt();
		if (numberOfGradDirection > 6)
		{
			tensorComputationPossible = true;
		}
		else//DTI data: can't get the right numberOfGradDirection from Philips Dicom Tag, calculate it myself when dealing with dwi data
		{
			DicomHelper::UpdateGradDirectionNumber();
		}
	}

	if (meta->HasAttribute(DICOM_DIFF_B_FACTOR))
	{
		//BvalueList.push_back(0);//Bug fix here, notify Wenxing later
		for (vtkDICOMDataElementIterator iter = meta->Begin(); iter != meta->End(); ++iter)
		{
			vtkDICOMTag tag = iter->GetTag();
			if ((tag == DICOM_DIFF_B_FACTOR) && (iter->IsPerInstance()))
			{
				//replace with below to fix possible multiple b value multiple direction error
				for (int i = 0; i < numberOfComponents; i = i + numberOfGradDirection)
				{
					BvalueList.push_back(meta->GetAttributeValue(i, DICOM_DIFF_B_FACTOR).AsFloat());
					std::cout << "BvalueList " << i / numberOfGradDirection << ": " << BvalueList.at(i / numberOfGradDirection) << std::endl;
				}			
				
				/*for (int i = 1; i < numberOfBValue; i++)
				{
					float bvalue = meta->GetAttributeValue(i, DICOM_DIFF_B_FACTOR).AsFloat();
					BvalueList.push_back(bvalue);
				}*/
			}
		}
	}

	if (tensorComputationPossible)
	{

		if (meta->HasAttribute(DICOM_IMAGE_ORIENTATION))
		{
			//std::string image_orientation;
			image_orientation = meta->GetAttributeValue(DICOM_IMAGE_ORIENTATION).GetCharData();
		}

	if (meta->HasAttribute(DICOM_IMGANG_RL) && meta->HasAttribute(DICOM_IMGANG_FH) && meta->HasAttribute(DICOM_IMGANG_AP))
	{
		ang[0] = meta->GetAttributeValue(DICOM_IMGANG_RL).AsDouble();
		ang[1] = meta->GetAttributeValue(DICOM_IMGANG_AP).AsDouble();
		ang[2] = meta->GetAttributeValue(DICOM_IMGANG_FH).AsDouble();
	}
	GetSliceToPatMatrix();

	if (meta->HasAttribute(DICOM_DIFF_DIRECT_AP) && meta->HasAttribute(DICOM_DIFF_DIRECT_FH) && meta->HasAttribute(DICOM_DIFF_DIRECT_RL))
	{
		gradientDirection = GradientDirectionContainerType::New();
		gradientDirection->Reserve(numberOfGradDirection);
		for (vtkDICOMDataElementIterator iter = meta->Begin(); iter != meta->End(); ++iter)
		{
			vtkDICOMTag tag = iter->GetTag();
			if ((tag == DICOM_DIFF_DIRECT_AP) && (iter->IsPerInstance()))
			{
				
				GradientDirectionContainerType::Iterator directionIterator = gradientDirection->Begin();
				for (int i = 1; i <= numberOfGradDirection; i++)
				{
					float AP_direction = meta->GetAttributeValue(i, DICOM_DIFF_DIRECT_AP).AsFloat();
					float FH_direction = meta->GetAttributeValue(i, DICOM_DIFF_DIRECT_FH).AsFloat();
					float RL_direction = meta->GetAttributeValue(i, DICOM_DIFF_DIRECT_RL).AsFloat();
					GradientDirectionType direction;
					direction[0] = RL_direction;
					direction[1] = AP_direction;
					direction[2] = FH_direction;
									
					float sum = abs(direction[0]) +abs(direction[1]) + abs(direction[2]);
					if (sum  < exp(-100))
					{
						IsoImageLabel = i; //the index in the component vector
					}
					else
					{
						directionIterator->Value() = direction;
						++directionIterator;
					}

					}
				}
			}
			if (IsoImageLabel > -1)
				numberOfGradDirection = numberOfGradDirection - 1;
			CalculateFinalHMatrix();
		}
	}
};

void DicomHelper::CalculateFinalHMatrix()
{
	finalH.set_size(numberOfGradDirection, 6);
	GradientDirectionContainerType::Iterator directionIterator = gradientDirection->Begin();
	for (int row = 0; row < numberOfGradDirection; row++)
	{
		GradientDirectionType direction;
		direction = directionIterator->Value();
		vnl_vector<double> newDirection(3);
		vnl_matrix<double> slice2PatInverse = vnl_matrix_inverse<double>(slice2PatMatrix);
		
		for (int i = 0; i < 3; i++)
		{
			double temp;
			temp = direction[0] * slice2PatInverse(i, 0) + direction[1] * slice2PatInverse(i, 1) + direction[2] * slice2PatInverse(i, 2);
			if (abs(temp) < exp(-8))
				newDirection(i) = 0.0;
			else
				newDirection(i) = temp;
			
		}
		finalH(row, 0) = newDirection[0] * newDirection[0];
		finalH(row, 1) = newDirection[1] * newDirection[1];
		finalH(row, 2) = newDirection[2] * newDirection[2];
		finalH(row, 3) = 2 * newDirection[0] * newDirection[1];
		finalH(row, 4) = 2 * newDirection[0] * newDirection[2];
		finalH(row, 5) = 2 * newDirection[1] * newDirection[2];
		directionIterator++;

		//cout << "newdirection: " << newDirection << endl;
	}
};

void DicomHelper::GetSliceToPatMatrix()
{
	vnl_matrix<double> slice_apat(3,3);
	if (strcmp(image_orientation, "TRANSVERSAL"))
	{
		slice_apat(0, 0) = 0.0; slice_apat(0, 1) = -1.0; slice_apat(0, 2) = 0.0;
		slice_apat(1, 0) = -1.0; slice_apat(1, 1) = 0.0; slice_apat(1, 2) = 0.0;
		slice_apat(2, 0) = 0.0; slice_apat(2, 1) =  0.0; slice_apat(2, 2) = 1.0;
		//cout << " transeveral " << endl;
	}
	else if (strcmp(image_orientation, "CORONAL"))
	{
		slice_apat(0, 0) = 0.0; slice_apat(0, 1) = 0.0; slice_apat(0, 2) = 1.0;
		slice_apat(1, 0) = -1.0; slice_apat(1, 1) = 0.0; slice_apat(1, 2) = 0.0;
		slice_apat(2, 0) = 0.0; slice_apat(2, 1) = 1.0; slice_apat(2, 2) = 0.0;
	}
	else if (strcmp(image_orientation, "SAGITTAL"))
	{
		slice_apat(0, 0) = 0.0; slice_apat(0, 1) = 0.0; slice_apat(0, 2) = 1.0;
		slice_apat(1, 0) = 0.0; slice_apat(1, 1) = -1.0; slice_apat(1, 2) = 0.0;
		slice_apat(2, 0) = -1.0; slice_apat(2, 1) = 0.0; slice_apat(2, 2) = 0.0;
	}

	vnl_matrix<double> apat_pat(3, 3);
	double	sx;  double	sy;  double	sz;
	double	cx;  double	cy;  double	cz;

	//if (right_handed)
	//{
	//	SGMAT_invert_vec(ang, &ang);
	//}

	sx = sin(-ang[0] * RAD); sy = sin(-ang[1] * RAD); sz = sin(-ang[2] * RAD);
	cx = cos(-ang[0] * RAD); cy = cos(-ang[1] * RAD); cz = cos(-ang[2] * RAD);
	
	apat_pat(0,0) = cy * cz;
	apat_pat(0,1) = -sz * cx + sx * sy * cz;
	apat_pat(0,2) = sx * sz + sy * cx * cz;

	apat_pat(1, 0) = sz * cy;
	apat_pat(1, 1) = cx * cz + sx * sy * sz;
	apat_pat(1, 2) = -sx * cz + sy * sz * cx;

	apat_pat(2, 0) = -sy;
	apat_pat(2, 1) = sx * cy;
	apat_pat(2, 2) = cx * cy;
	//cout << "apat to pat matrix " << endl;
	//cout << apat_pat << endl;

	slice2PatMatrix.set_size(3, 3);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3;j++)
		{
			slice2PatMatrix(i, j) = apat_pat(0, j)*slice_apat(i, 0) + apat_pat(1, j)*slice_apat(i, 1) + apat_pat(2, j)*slice_apat(i, 2);
		}
		
	}
	vnl_matrix<double> multiplier(3, 3);
	multiplier(0, 0) = 0.0; multiplier(0, 1) = -1.0; multiplier(0, 2) = 0.0;
	multiplier(1, 0) = -1.0; multiplier(1, 1) = 0.0; multiplier(1, 2) = 0.0;
	multiplier(2, 0) = 0.0; multiplier(2, 1) = 0.0; multiplier(2, 2) = 1.0;

	slice2PatMatrix = slice2PatMatrix.transpose()*multiplier;
};


void DicomHelper::UpdateGradDirectionNumber()
{
	numberOfGradDirection = numberOfBValue > 1 ? (numberOfComponents - 1) / (numberOfBValue - 1) : numberOfGradDirection;
};

DicomHelper::DicomHelper(vtkStringArray* Files)
{
	//cout << "direcotry name" << DirectoryName << endl;
	//dDir = vtkSmartPointer<vtkDICOMDirectory>::New();
	//dDir->SetDirectoryName(DirectoryName);
	//dDir->SetScanDepth(1);
	//dDir->Update();	
	//FileNamesForDiffusionSeries = Files;
	DicomReader = vtkSmartPointer<vtkDICOMReader>::New();
	DicomReader->SetFileNames(Files);

	DicomHelper::DicomInfo();
	//DicomReader->Update();
	//numberOfComponents = DicomReader->GetOutput()->GetNumberOfScalarComponents();
	//DicomReader->GetOutput()->GetDimensions(imageDimensions);

};
