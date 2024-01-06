from PyQt6.QtCore import pyqtSlot, QSettings, QTimer, QUrl, QDir, Qt
from PyQt6.QtWidgets import QMainWindow, QFileDialog, QMessageBox, QApplication, QDockWidget, QPlainTextEdit
from PyQt6.QtWebEngineWidgets import QWebEngineView

import sys
import subprocess
import signal
import logging
import threading

global logger, app, ui
    
#define UI elements
class LoggerDock(QDockWidget):

    def __init__(self, *args):
        super(LoggerDock, self).__init__(*args)
        self.textview = QPlainTextEdit(self)
        self.textview.setReadOnly(True)
        self.setWidget(self.textview)

    @pyqtSlot(str)
    def log(self, message):
        self.textview.appendPlainText(message)
        
class mainGUI(QMainWindow):

    def __init__(self, parent=None, homepage=None, *args, **kwargs):
        super(mainGUI, self).__init__(*args, **kwargs)
        self.homepage = homepage
        self.windows = []

        self.loggerdock = LoggerDock("Log Message", self)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.loggerdock)
        
        settings = QSettings()
        val = settings.value("net.fishandwhistle/JupyterQt/geometry", None)
        if val is not None:
            self.restoreGeometry(val)

        self.basewebview = JupyterView(self, main=True)
        self.setCentralWidget(self.basewebview)
        #self.windows.append(self.basewebview)
        #self.tabs = QTabWidget(self)
        #self.tabs.setTabsClosable(True)
        #self.tabs.setMovable(True)
        #self.tabs.tabCloseRequested.connect(self.destroyBrowserTab)
        
        #self.setCentralWidget(self.tabs)
        self.show()
        
        QTimer.singleShot(0, self.initialload)

    @pyqtSlot()
    def initialload(self):
        if self.homepage:
            self.basewebview.load(QUrl(self.homepage))
        self.show()

    def closeEvent(self, event):
        if self.windows:
            for i in reversed(range(len(self.windows))):
                w = self.windows.pop(i)
                w.close()
            event.accept()
        else:
            event.accept()

        #save geometry
        settings = QSettings()
        settings.setValue("net.fishandwhistle/JupyterQt/geometry", self.saveGeometry())

class JupyterView(QWebEngineView):

    def __init__(self, mainwindow, main=False):
        super(JupyterView, self).__init__(None)
        self.parent = mainwindow
        self.main = main
        self.loadedPage = None

    @pyqtSlot(bool)
    def onpagechange(self, ok):
        log("on page change: %s, %s" % (self.url(), ok))
        if self.loadedPage is not None:
            log("disconnecting on close signal")
            self.loadedPage.windowCloseRequested.disconnect(self.close)
        self.loadedPage = self.page()
        log("connecting on close signal")
        self.loadedPage.windowCloseRequested.connect(self.close)

    def createWindow(self, windowtype):
        v = JupyterView(self.parent)
        windows = self.parent.windows
        windows.append(v)
        v.show()
        return v

    def closeEvent(self, event):
        if self.loadedPage is not None:
            log("disconnecting on close signal")
            self.loadedPage.windowCloseRequested.disconnect(self.close)

        if not self.main:
            if self in self.parent.windows:
                self.parent.windows.remove(self)
            log("Window count: %s" % (len(self.parent.windows)+1))
        event.accept()


def initializeJupyterLab():
    #start jupyter notebook and wait for line with the web address
    log("Starting Jupyter notebook process")
    notebookp = startnotebook()

    log("Waiting for server to start...")
    webaddr = None
    while webaddr is None:
        line = str(notebookp.stderr.readline())
        log(line)
        if "http://" in line:
            start = line.find("http://")
            end = line.find("\n", start+len("http://"))
            webaddr = line[start:end-2]
            print(webaddr)
        log("Server found at %s, migrating monitoring to listener thread" % webaddr)

    #pass monitoring over to child thread
    def process_thread_pipe(process):
        while process.poll() is None: #while process is still alive
            log(str(process.stderr.readline()))

    notebookmonitor = threading.Thread(name="Notebook Monitor", target=process_thread_pipe,
                                       args = (notebookp,))
    notebookmonitor.start()
    return notebookp, webaddr
    

def startnotebook(notebook_executable="jupyter-lab", port=8888, directory='/Users/mb/Documents/Syntuitio/AudreyHay/PlanB/python/'): #QDir.homePath()
    return subprocess.Popen([notebook_executable,
                            "--port=%s" % port, "--browser=n", "-y",
                            "--notebook-dir=%s" % directory], bufsize=1,
                            stderr=subprocess.PIPE)


def handle_exception(exc_type, exc_value, exc_traceback):
    """Handle uncaught exceptions and print in logger."""
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    
def log(message):
    logger.debug(message)
    try:
        ui.loggerdock.log(message)
    except NameError:
        pass
    
if __name__ == '__main__':
    global logger, app, ui#, pr
    
    # logging
    if getattr(sys, 'frozen', False):
        fnLog = os.path.join(sys._MEIPASS, 'interfaceJupyter.log')
    else:
        fnLog = 'interfaceJupyter.log'
            
    logger = logging.getLogger("my-app")
    logger.setLevel(logging.DEBUG)
    
    fh = logging.FileHandler(fnLog)
    fh.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        '%(asctime)s :: %(filename)s :: %(funcName)s :: line %(lineno)d ::'
        ' %(levelname)s :: %(message)s',
        datefmt='%Y-%m-%d:%H:%M:%S')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    sys.excepthook = handle_exception
    
    notebookp, webaddr=initializeJupyterLab()

    #setup application
    log("Setting up GUI")
    app = QApplication(sys.argv)
    app.setApplicationName("JupyterQt")
    app.setOrganizationDomain("fishandwhistle.net")

    #setup webview
    ui = mainGUI(None, homepage=webaddr)

    log("Starting Qt Event Loop")
    result = app.exec()

    log("Sending interrupt signal to jupyter-notebook")
    notebookp.send_signal(signal.SIGINT)
    try:
        log("Waiting for jupyter to exit...")
        notebookp.wait(10)
    except subprocess.TimeoutExpired:
        log("control c timed out, killing")
        notebookp.kill()

    log("Exited.")
    sys.exit(result)
