#include "MantidAPI/NotebookWriter.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/MantidVersion.h"

namespace Mantid {
namespace API {

namespace {
    Mantid::Kernel::Logger g_log("NotebookWriter");
}

NotebookWriter::NotebookWriter() : m_cell_buffer(Json::arrayValue)
{
  // Add header comments and code
  headerComment();
  headerCode();
}

  /**
 * Add a code cell to the buffer of cells to write to the notebook
 *
 * @param array_code :: Json array of strings containing python code
   * for the code cell
 */
void NotebookWriter::codeCell(Json::Value array_code) {

  Json::Value cell_data;
  const Json::Value empty = Json::Value(Json::ValueType::objectValue);

  cell_data["cell_type"] = "code";
  cell_data["collapsed"] = false;
  cell_data["input"] = array_code;
  cell_data["language"] = "python";
  cell_data["metadata"] = empty;
  cell_data["outputs"] = Json::Value(Json::arrayValue);

  m_cell_buffer.append(cell_data);
}

  /**
 * Add a code cell to the buffer of cells to write to the notebook
 *
 * @param string_code :: string containing the python for the code cell
 */
std::string NotebookWriter::codeCell(std::string string_code) {

  Json::Value cell_data;
  const Json::Value empty = Json::Value(Json::ValueType::objectValue);

  cell_data["cell_type"] = "code";
  cell_data["collapsed"] = false;
  cell_data["input"] = string_code;
  cell_data["language"] = "python";
  cell_data["metadata"] = empty;
  cell_data["outputs"] = Json::Value(Json::arrayValue);

  m_cell_buffer.append(cell_data);
  Json::StyledWriter writer;
  return writer.write(cell_data);
}

  /**
* Add a markdown cell to the buffer of cells to write to the notebook
*
* @param string_array :: json array of strings containing the python
   * code for the code cell
*/
void NotebookWriter::markdownCell(Json::Value string_array) {

  Json::Value cell_data;
  const Json::Value empty = Json::Value(Json::ValueType::objectValue);

  cell_data["cell_type"] = "markdown";
  cell_data["metadata"] = empty;
  cell_data["source"] = string_array;

  m_cell_buffer.append(cell_data);
}

  /**
* Add a markdown cell to the buffer of cells to write to the notebook
*
* @param string_text :: string containing the python code for the code cell
*/
std::string NotebookWriter::markdownCell(std::string string_text) {

  Json::Value cell_data;
  const Json::Value empty = Json::Value(Json::ValueType::objectValue);

  cell_data["cell_type"] = "markdown";
  cell_data["metadata"] = empty;
  cell_data["source"] = string_text;

  m_cell_buffer.append(cell_data);
  Json::StyledWriter writer;
  return writer.write(cell_data);
}

  /**
* Add a markdown cell of information for the user to the buffer of cells to
   * write to the notebook
*/
void NotebookWriter::headerComment() {

  Json::Value strings(Json::arrayValue);
  strings.append(Json::Value("This IPython Notebook was automatically generated by MantidPlot, version: "));
  strings.append(Json::Value(Mantid::Kernel::MantidVersion::version()));
  strings.append(Json::Value("\n"));
  strings.append(Json::Value(Mantid::Kernel::MantidVersion::releaseNotes()));
  strings.append(Json::Value("\n\nThe following information may be useful:\n"
                                 "* [Mantid Framework Python API Reference]"
                                 "(http://docs.mantidproject.org/nightly/api/python/index.html)\n"
                                 "* [IPython Notebook Documentation](http://ipython.org/ipython-doc/stable/notebook/)\n"
                                 "* [matplotlib Documentation](http://matplotlib.org/contents.html)\n\n"
                                 "Help requests and bug reports should be submitted to the [Mantid forum.]"
                                 "(http://forum.mantidproject.org)"));

  markdownCell(strings);
}

  /**
* Add code cells to the buffer of cells to write to the notebook
   * These are to import Mantid and matplotlib, and to warn the
   * user if the version of Mantid being used does not match the
   * version which generated the notebook.
*/
void NotebookWriter::headerCode() {

  Json::Value import_mantid(Json::arrayValue);

  import_mantid.append(Json::Value("import sys\n"
                                   "import os\n"
                                   "\n"
                                   "#Find where Mantid is installed and tell python.\n"
                                   "def find_path():\n"
                                   "    for r,d,f in os.walk(os.path.abspath(os.sep)):\n"
                                   "        for files in f:\n"
                                   "            if files == 'MantidPlot.exe' or files == 'MantidPlot' or files == 'MantidPlot.app':\n"
                                   "                return r\n"
                                   "\n"
                                   "mantidPath = 'C://MantidInstall/bin'\n"
                                   "if os.path.isdir(mantidPath): sys.path.append(mantidPath)\n"
                                   "else: sys.path.append(find_path())\n"
                                   "\n"
                                   "#We can now import Mantid's Python API\n"
                                   "from mantid.simpleapi import *\n"
                                   "#If import fails then replace path on line \"mantidPath = ...\" with correct path to the MantidPlot executable"));

  codeCell(import_mantid);

  Json::Value check_version(Json::arrayValue);

  check_version.append(Json::Value("import warnings"));
  check_version.append(Json::Value("import mantid.kernel"));
  check_version.append(Json::Value("# Check if the version of Mantid being used matches"
                                       " the version which created the notebook."));
  check_version.append(Json::Value(std::string("if \"") + Mantid::Kernel::MantidVersion::version() +
                                   "\" != mantid.kernel.version_str(): warnings.warn(\"Version of Mantid"
                                       " being used does not match version which created the notebook.\")"));

  codeCell(check_version);

  Json::Value import_matplotlib(Json::arrayValue);

  import_matplotlib.append(Json::Value("#Import matplotlib's pyplot interface under the name 'plt'\n"));
  import_matplotlib.append(Json::Value("import matplotlib.pyplot as plt\n\n"));

  import_matplotlib.append(Json::Value("#Some magic to tell matplotlib how to behave in IPython Notebook. "
                                         "Use '%matplotlib nbagg' for interactive plots, if available.\n"));
  import_matplotlib.append(Json::Value("%matplotlib inline"));

  codeCell(import_matplotlib);
}

  /**
* Create a Json value containing the whole notebook
 * @return a Json value containing the whole notebook
*/
Json::Value NotebookWriter::buildNotebook() {

  Json::Value output;
  const Json::Value empty = Json::Value(Json::ValueType::objectValue);

  Json::Value worksheet;
  worksheet["cells"] = m_cell_buffer;
  worksheet["metadata"] = empty;

  Json::Value worksheet_arr(Json::arrayValue);
  worksheet_arr.append(worksheet);

  Json::Value meta_name;
  meta_name["name"] = "Mantid Notebook";
  output["metadata"] = meta_name;
  output["nbformat"] = 3;
  output["nbformat_minor"] = 0;
  output["worksheets"] = worksheet_arr;

  return output;
}

  /**
* Create a formatted string of Json which describes a notebook
   * @return a formatted string of the Json which describes
   * the whole notebook
*/
std::string NotebookWriter::writeNotebook() {

  const Json::Value root = buildNotebook();

  Json::StyledWriter writer;
  std::string output_string = writer.write(root);

  return output_string;
}

}
}
