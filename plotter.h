#ifndef PLOTTER_H
#define PLOTTER_H

#include <QMap>
#include <QPixmap>
#include <QVector>
#include <QWidget>
#include <QTimer>
#include <QMutex>
#include <QDateTime>
#include <QFile>
#include "DAQReader.h"
#include "QuadFilter.h"

class QToolButton;
class PlotSettings;

class Plotter : public QWidget
{
    Q_OBJECT

    public:
        Plotter(QWidget *parent = 0);

        QSize minimumSizeHint() const;
        QSize sizeHint() const;

    public slots:
        void zoomIn();
        void zoomOut();
        void toggleRecording();
        void newData();
        void daqError(const QString& errorMessage);
        void startedRecording();
        void stoppedRecording();
        void newDocument();
        void open();
        void save();
        void settings();
        void filtering();

    protected:
        void paintEvent(QPaintEvent *event);
        void resizeEvent(QResizeEvent *event);
        void mousePressEvent(QMouseEvent *event);
        void mouseMoveEvent(QMouseEvent *event);
        void mouseReleaseEvent(QMouseEvent *event);
        void keyPressEvent(QKeyEvent *event);
        void wheelEvent(QWheelEvent *event);
        void closeEvent(QCloseEvent* event);

    private:
        bool offerToSave();
        void clearPlot();
        void updateRubberBandRegion();
        void updateFilteredData();
        void refreshPixmap();
        void drawGrid(QPainter *painter);
        void drawCurves(QPainter *painter);
        void updateSettings();

        enum { Margin = 50 };

        DAQSettings daqSettings;

        double fileSamplingRate;

        QToolButton *newButton;
        QToolButton *openButton;
        QToolButton *saveButton;
        QToolButton *settingsButton;
        QToolButton *filteringButton;
        QToolButton *recordButton;
        QToolButton *zoomInButton;
        QToolButton *zoomOutButton;
        typedef QMap<int, QVector<QPointF> > CurveMap;
        CurveMap curveMap;
        typedef QMap<int, QVector<QuadFilter> > FilterMap;
        FilterMap filterMap;
        CurveMap filteredCurveMap;
        QVector<PlotSettings> zoomStack;
        int curZoom;
        bool rubberBandIsShown;
        QRect rubberBandRect;
        QPixmap pixmap; // offscreen rendering buffer
        bool saved;
        QDateTime startTime;
        DAQReader daqReader;
        QString filename;
        double traceOffset;
        QFile sharedTimestamp;
        uchar* sharedTimestampMemMap;

#ifdef Q_WS_MAC
        bool recording;
#endif
};

class PlotSettings
{
    public:
        PlotSettings();

        void scroll(double dx, double dy);
        void adjust();
        double spanX() const { return maxX - minX; }
        double spanY() const { return maxY - minY; }

        double minX;
        double maxX;
        int numXTicks;
        double minY;
        double maxY;
        int numYTicks;

        bool includesRightEdge;

    private:
        static void adjustAxis(double &min, double &max, int &numTicks);
};

#endif
