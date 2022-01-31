# Going from Qt Designer to code
- Save layout from Qt Designer to a `.ui` file
- Run `pyside2-uic` on the `.ui` file and direct the output to a python file. This generates a class
  containing all the UI elements.
  * For example, if the UI elements are saved in `mockup.ui`, then we can do:
  ```
  pyside2-uic mockup.ui > ui_mainwindow.py
  ```
- Import the UI class into the main python script, e.g.: `from ui_mainwindow import Ui_MainWindow`
- You can now use the widgets in the layout via the following pattern:

  ```
  class my_class(...):
    self.ui = Ui_MainWindow()
    self.ui.setupUi(self)
    ...
    # Register a function with the "start_button" widget
    self.ui.start_button.clicked.connect(some_function)
  ```
