/*===================================================================

DWI Core UI

===================================================================*/

#ifndef DiffusionCore_h
#define DiffusionCore_h

//// include ctk
//#include <ctkDICOMDatabase.h>
//#include <ctkDICOMIndexer.h>

//include QT
#include <QHash>
#include <QLabel>
#include <QProgressDialog>
#include <QString>
#include <QStringList>
#include <QVariant>
#include <QWidget>
#include <qdebug.h>

class ctkFileDialog;

namespace Ui{
	class DiffusionModule;
}
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

	///// @brief
	//void OnViewButtonClicked();

	///// @brief
	//void OnStartDicomImport(const QString&);

	///// @brief
	///// In this slot the models database with new imports is set.
	///// This causes a model update.
	//void OnFinishedImport();

	//void OnSeriesSelectionChanged(const QStringList &s);

	protected slots:

	/// \brief 
	//void OnScanDirectory();

protected:

	/// \brief 
	
	//QStringList GetFileNamesFromIndex();

	/// \brief 
	
	//void SetupImportDialog();

	/// \brief 
	
	//void SetupProgressDialog(QWidget* parent);

	//ctkDICOMDatabase* m_ExternalDatabase;
	//ctkDICOMIndexer* m_ExternalIndexer;
	//ctkFileDialog* m_ImportDialog;

	//QLabel* m_ProgressDialogLabel;
	//QProgressDialog* m_ProgressDialog;
	//QString m_LastImportDirectory;

	Ui::DiffusionModule* m_Controls;

};



#endif //

