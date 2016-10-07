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

	///// @brief emitted when import into database is finished.
	//void SignalStartDicomImport(const QStringList&);


	///// @brief emitted when view button is clicked.
	//void SignalDicomToDataManager(QHash<QString, QVariant>);



	public slots:

	/// @brief Called when dicom files were input.
	void OnImageFilesLoaded(const QStringList&);

	protected slots:

	///// @brief
	///// In this slot, implement the adc calculation
	///// This causes a model update.
	void calcADC(bool toggle);

	///// @brief
	///// In this slot, implement the cDWI calculation and image rendering
	///// 
	void calcCDWI(bool toggle);

	void calcEADC(bool toggle);

	///// @brief
	///// In this slot, implement the cdwi calculation with input bvalue
	///// 
	void cDWI(double bvalue);

	///// @brief
	///// In this slot, implement the filtering with input threshhold
	///// 
	void adc(double threshhold);

protected:

	Ui::DiffusionModule* m_Controls;
	DicomHelper *m_DicomHelper;//initialization? 
	//QString m_LastImportDirectory;
	void DisplayDicomInfo(vtkSmartPointer <vtkImageData> imageData);
	void SourceImageViewer2D(vtkSmartPointer <vtkImageData>, QVTKWidget *qvtkWidget);
	void QuantitativeImageViewer2D(vtkSmartPointer <vtkImageData>, QVTKWidget *qvtkWidget);
	void TestCallbackFunc(vtkObject *caller, long unsigned int eventId, void *clientData, void* callData);


protected:
	vtkSmartPointer< vtkImageData > sourceImage;
	vtkSmartPointer <vtkRenderWindowInteractor> m_RenderWindowInteractor;
	vtkSmartPointer< vtkImageData > cacheImage;
	int m_SourceImageCurrentSlice;
	int m_QuantitativeImageCurrentSlice;
	double m_MaskThreshold;
	double m_ComputedBValue;
};


//class DiffusionCore : public QVTKWidget
//{
//	protected:
//}

#endif //

