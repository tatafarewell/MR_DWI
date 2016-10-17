/*===================================================================

DWI Core UI

===================================================================*/

#ifndef DiffusionCore_h
#define DiffusionCore_h

// Headers in this module
// #include <DicomHelper.h>

//include QT
//#include <QHash>
//#include <QLabel>
//#include <QProgressDialog>
//#include <QString>
//#include <QStringList>
//#include <QVariant>

// include VTK
#include <vtkImageData.h>
#include <vtkRenderWindowInteractor.h>
#include "vtkSmartPointer.h"

#include <QWidget>
#include <qdebug.h>

//#include "QVTKWidget.h"

class ctkFileDialog;

namespace Ui{
	class DiffusionModule;
}
class DicomHelper;
class QVTKWidget;
/**
* \brief DiffusionCore is a QWidget providing functionality for diffusion weighted image calculation.
*
* \sa 
* \ingroup 
*/
#define MAXWINDOWNUMBER 5 //difine the busiest viewing window layout as 5 by 5
class DiffusionCore : public QWidget
{
	// this is needed for all Qt objects that should have a Qt meta-object
	// (everything that derives from QObject and wants to have signal/slots)
	Q_OBJECT

public:

	//static const std::string Widget_ID;

	/**
	* \brief DiffusionCore(QWidget *parent) constructor.
	*
	* \param parent is a pointer to the parent widget
	*/
	DiffusionCore(QWidget *parent);

	/**
	* \brief DiffusionCore destructor.
	*/
	virtual ~DiffusionCore();

	/**
	* \brief CreateQtPartControl(QWidget *parent) sets the view objects from ui_DiffusionModule.h.
	*
	* \param parent is a pointer to the parent widget
	*/
	virtual void CreateQtPartControl(QWidget *parent);

	/**
	* \brief Initializes the widget. This method has to be called before widget can start.
	*/
	void Initialize();


signals:

	///// @brief emitted when dicomdata is imported.
	void SignalDicomLoaded(bool);

	public slots:

	///// @brief 
	///// In this slot,  render input dicom files.
	///// 
	void OnImageFilesLoaded(const QStringList&);

	protected slots:

	///// @brief
	///// In this slot, call adc calculation and render the image
	///// 
	void onCalcADC(bool toggle);

	///// @brief
	///// In this slot, call cdwi calculation and render the image
	///// 
	void onCalcCDWI(bool toggle);

	///// @brief
	///// In this slot, call eADC calculation and render the image
	///// 
	void onCalcEADC(bool toggle);

	///// @brief
	/////In this slot, call FA calculation and render the image
	///// 
	void onCalcFA(bool toggle);

	///// @brief
	///// In this slot, call colorFA calculation and render the image
	/////
	void onCalcColorFA(bool toggle);

	///// @brief
	///// In this slot, call IVIM calculation and render the image
	///// 
	void onCalcIVIM(bool toggle);

	///// @brief
	///// In this slot, update the cdwi image with input bvalue
	///// 
	void onBSlide(double bvalue);

	///// @brief
	///// In this slot, update the adc image with filtering using input threshhold
	///// 
	void onThreshSlide(double threshhold);

	//void outputMOUSEevent(QMouseEvent*);

protected:
	enum imageType
	{
		ORIGINAL = 0,
		ADC = 1,
		CDWI = 2,
		EADC = 3,
		FA = 4,
		CFA = 5,
		IVIM = 6,
		NOIMAGE = 100
	};

	Ui::DiffusionModule* m_Controls;
	DicomHelper *m_DicomHelper;//initialization? 
	//QString m_LastImportDirectory;

	void DisplayDicomInfo(vtkSmartPointer <vtkImageData> imageData);
	void SourceImageViewer2D(vtkSmartPointer <vtkImageData>, QVTKWidget *qvtkWidget);
	void QuantitativeImageViewer2D(vtkSmartPointer <vtkImageData>, QVTKWidget *qvtkWidget, std::string imageLabel);
	void TestCallbackFunc(vtkObject *caller, long unsigned int eventId, void *clientData, void* callData);
	void ShareWindowEvent();

	void AdcCalculator(vtkSmartPointer <vtkImageData> imageData);
    void FaCalculator(vtkSmartPointer <vtkImageData> imageData);
	void ColorFACalculator(vtkSmartPointer <vtkImageData> imageData);
	void EAdcCalculator(vtkSmartPointer <vtkImageData> imageData);	
	void CDWICalculator(vtkSmartPointer <vtkImageData> imageData);

	void ui_InsertWindow(int& rowInd, int& colInd, QVTKWidget *vtkWindow, imageType imageLabel);
	void ui_RemoveWindow(imageType a);
	bool ui_IsWdWSquare();
	void ui_dumpWindow(int row, int col);
	void ui_findWdw(imageType imageLabel, int& row, int& col);
	

protected:
	vtkSmartPointer < vtkImageData > sourceImage;
	vtkSmartPointer <vtkRenderWindowInteractor> m_RenderWindowInteractor;
	vtkSmartPointer< vtkImageData > cacheImage;
	int m_SourceImageCurrentSlice;
	int m_QuantitativeImageCurrentSlice;
	//bool m_SourceImageIntera
	double m_MaskThreshold;
	double m_ComputedBValue;	
	std::vector< std::vector<int> > layoutTable; // Table used for tracing window content

};



#endif //

