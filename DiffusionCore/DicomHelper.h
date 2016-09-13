#include <vtkDICOMDirectory.h>
#include <vtkSmartPointer.h>
#include <vtkDICOMReader.h>
#include <vnl/vnl_matrix.h>

#include <DTISolution.h>

#include <iostream>
#include <vector>

#define PI 3.141592653589793
#define RAD  (PI/180.0)


#define DICOM_DIFF_B_FACTOR       vtkDICOMTag(0x2001, 0x1003)
#define DICOM_DIFF_DIRECT         vtkDICOMTag(0x2001, 0x1004)
#define DICOM_SCALE_SLOP          vtkDICOMTag(0x2005, 0x100E)
#define DICOM_SCALE_INTERPCEPT    vtkDICOMTag(0x2005, 0x100D)
#define DICOM_NO_DIFF_BVAL        vtkDICOMTag(0x2005, 0x1414)
#define DICOM_NO_DIFF_GRAD_ORIENT vtkDICOMTag(0x2005, 0x1415)
#define DICOM_DIFF_DIRECT_RL      vtkDICOMTag(0x2005, 0x10b0) //x direction
#define DICOM_DIFF_DIRECT_AP      vtkDICOMTag(0x2005, 0x10b1) //z direction
#define DICOM_DIFF_DIRECT_FH      vtkDICOMTag(0x2005, 0x10b2) //y direction
#define DICOM_PLANE_ORIENTATION   vtkDICOMTag(0x0020, 0x0037) //double
#define DICOM_IMGANG_RL           vtkDICOMTag(0x2005, 0x1002)//RL
#define DICOM_IMGANG_FH           vtkDICOMTag(0x2005, 0x1001)//FH
#define DICOM_IMGANG_AP           vtkDICOMTag(0x2005, 0x1000)//AP
#define DICOM_IMAGE_ORIENTATION   vtkDICOMTag(0x2001, 0x100B) //trans cor sag



class DicomHelper
{
public:
	DicomHelper(char *DirecotyName);
	void Dicomread();
	void DicomHelper::DicomInfo();
	vtkSmartPointer<vtkDICOMReader> DicomReader;
	~DicomHelper() {};
	void DicomHelper::RemoveIsoImage();

	float scaleSlope;
	float scaleIntercept;
	int numberOfGradDirection;
	int numberOfBValue;
	int IsoImageLabel = -1; // -1 means no isotropic image

	std::vector <float> BvalueList;
	std::vector <float> directionLabel;
	vnl_matrix<float> finalH;

	const char *image_orientation;
	itk::Vector<double, 3> ang;
	vnl_matrix<double> slice2PatMatrix;

	GradientDirectionContainerType::Pointer gradientDirection;

private:
	int DiffusionSeryNumber;
	vtkSmartPointer<vtkDICOMDirectory> dDir;

	vtkStringArray *FileNamesForDiffusionSeries;
	int GetDiffusionDataset(char *DirectoryName);

	void DicomHelper::CalculateFinalHMatrix();
	void DicomHelper::GetSliceToPatMatrix();


};