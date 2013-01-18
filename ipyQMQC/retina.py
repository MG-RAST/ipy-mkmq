#!/usr/bin/env python

import IPython.core.display
import json
import ipyTools

class Retina(object):
    def __init__(self, action='none', debug=False):
        self.action = action
        self.debug  = debug
        self.rjs = ipyTools.Ipy.RETINA_URL+'js/'
        self.rcss   = [ ipyTools.Ipy.RETINA_URL+'css/bootstrap.min.css' ]
        self.rlibs  = [ self.rjs+'bootstrap.min.js',
                        self.rjs+'retina.js',
                        self.rjs+'stm.js',
                        self.rjs+'ipy.js' ]
        self.renderer_resource = ipyTools.Ipy.RETINA_URL+"renderers/";
        src = """
			(function(){
				Retina.init( { library_resource: '"""+self.rjs+"""'});
			})();
		"""
        IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src, lib=self.rlibs, css=self.rcss))
    
    def graph(self, width=800, height=400, btype="column", target="", data=None, title="", x_labels=[], x_title="", y_title="", show_legend=False, legend_position='left', title_color="black", x_title_color="black", y_title_color="black", x_labels_rotation="0", x_tick_interval=0, y_tick_interval=30, x_labeled_tick_interval=1, y_labeled_tick_interval=5, default_line_color="black", default_line_width=1, chartArea=None, legendArea=None, onclick=None):
        """Graph Renderer
  
  Displays a graph of pie / bar charts with an optional legend.
  
  Options
  
  btype (STRING)
      Defines the display type of the graph, can be one of
        pie
        column
        stackedColumn
        row
        stackedRow
        line
      Default is column.
  
  title (STRING)
      Title string written at the top of the graph
  
  title_color (CSS Color Value)
      Color of the title text. Default is black.
  
  x_title (STRING)
      Title written below the x-axis.
  
  y_title (STRING)
      Title written to the left of the y-axis.
  
  x_title_color (CSS Color Value)
      Color of the x-axis title string. Default is black.
  
  y_title_color (CSS Color Value)
      Color of the y-axis title string. Default is black.
  
  x_labels (ARRAY of STRING)
      List of the labels at the ticks of the x-axis.
  
  x_labels_rotation (STRING)
      A string representing the number of degrees to rotate the labels on the x-axis. Default is 0.
  
  y_labels (ARRAY of STRING)
      List of the labels at the ticks of the y-axis. If no list is passed will use the y-valus.
  
  x_tick_interval (INT)
      Determines how many ticks are actually drawn on the x-axis. Default is 0.
  
  y_tick_interval (INT)
      Determines how many ticks are actually drawn on the y-axis. Default is 30.
  
  x_labeled_tick_interval (INT)
      Determines which ticks on the x-axis get labels. Default is 1.
  
  y_labeled_tick_interval (INT)
      The number of y-axis ticks that get labels. Default is 5.
  
  default_line_color (CSS Color Value)
      Determines the color of lines if not specified for an individual line. Default is black.
  
  default_line_width (INT)
      Number of pixels lines should be wide if not specified for an individual line. Default is 1.
  
  show_legend (BOOLEAN)
      Turns the display of the legend on / off. Default ist true.
  
  legend_position (STRING)
      Can be one of
        left
        right
        top
        bottom
  
  chartArea (ARRAY of FLOAT)
     The values passed correspond to the left, top, width and height margin of the chart area respectively. The position is relative to the top left corner of the containing div. Values less than 1 are interpreted as fractions. Values greater than 1 are interpreted as absolute pixel values.
  
  legendArea (ARRAY of FLOAT)
      If this parameter is set, the legend_position parameter will not be used. Instead pass an array of floats. The values correspond to the left, top, right and bottom margin of the legend area respectively. The position is relative to the top left corner of the containing div. Values less than 1 are interpreted as fractions. Values greater than 1 are interpreted as absolute pixel values.
  
  width (INT)
      The width of the graph in pixel (including legend).
  
  height (INT)
      The height of the graph in pixel (including legend).
  
  data (ARRAY of OBJECT)
      List of data series. Each series has a name and a data attribute. The data attribute is a list of y-values for the series.
  
  onclick (FUNCTION)
      The passed function will be called when a bar / pie slice is clicked. It will receive an object with the attributes
        series - the name of the series this bar belongs to
        value  - the value of the bar
        label  - the label of the bar
        item   - the svg element that was clicked
        index  - the zero based index of this bar within its series
        series_index - the zero based index of this series"""
        if not target:
            target = 'div_'+ipyTools.random_str()
        html = "<div id='%s'></div>"%(target)
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
        if len(x_labels) == 0:
            x_labels = [""]
        if data is None:
            title = "Browser Usage"
            x_labels = "['2005','2006','2007','2008']"
            data = "Retina.RendererInstances.graph[0].exampleData()"
            onclick = "'clickedCell = '+ JSON.stringify(params)"
        else:
            data = json.dumps(data)
            x_labels = json.dumps(x_labels)

        opt = "width: %d, height: %d, type: '%s', target: document.getElementById('%s'), data: %s, title: '%s', x_labels: %s, x_title: '%s', y_title: '%s', show_legend: %s, legend_position: '%s', title_color: '%s', x_title_color: '%s', y_title_color: '%s', x_labels_rotation: '%s', x_tick_interval: %f, y_tick_interval: %f, x_labeled_tick_interval: %f, y_labeled_tick_interval: %f, default_line_color: '%s', default_line_width: %d"%(width, height, btype, target, data, title, x_labels, x_title, y_title, self._bool(show_legend), legend_position, title_color, x_title_color, y_title_color, x_labels_rotation, x_tick_interval, y_tick_interval, x_labeled_tick_interval, y_labeled_tick_interval, default_line_color, default_line_width)
        
        if chartArea:
            opt += ", chartArea: "+json.dumps(chartArea)
        if legendArea:
            opt += ", legendArea: "+json.dumps(legendArea)
        if onclick:
            onclick = ", onclick: function(params){ipy.write_cell(ipy.add_cell(),"+onclick+");}"
        else:
            onclick = ""
        
        src = """
			(function(){
				Retina.add_renderer({"name": "graph", "resource": '""" + self.renderer_resource + """', "filename": "renderer.graph.js" });
				Retina.load_renderer("graph").then( function () { Retina.Renderer.create('graph', {""" + opt + onclick + """}).render(); });
                        })();
		"""
        if self.debug:
            print src
        else:
            IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))
    
    def plot(self, width=800, height=400, target="", data=None, title="", show_legend=True, legend_position='left'):
        if not target:
            target = 'div_'+ipyTools.random_str()
        html = "<div id='%s'></div>"%(target)
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
        if data is None:
            title = "Sine"
            data = "Retina.RendererInstances.plot[0].exampleData()"
        else:
            data = json.dumps(data)
        
        opt = "width: %d, height: %d, target: document.getElementById('%s'), data: %s, title: '%s', show_legend: %s, legend_position: '%s'"%(width, height, target, data, title, self._bool(show_legend), legend_position)
        src = """
			(function(){
				Retina.add_renderer({"name": "plot", "resource": '""" + self.renderer_resource + """', "filename": "renderer.plot.js" });
				Retina.load_renderer("plot").then( function () { Retina.Renderer.create('plot', {""" + opt + """}).render(); });
                        })();
		"""
        if self.debug:
            print src
        else:
            IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))
    
    def paragraph(self, width="span12", target="", data=None, title_color='black', header_color='black', text_color='black', raw=False):
        if not target:
            target = 'div_'+ipyTools.random_str()
        html = "<div id='%s'></div>"%(target)
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
        if data is None:
            data = "Retina.RendererInstances.paragraph[0].exampleData()"
        else:
            data = json.dumps(data)
        
        opt = "width: '%s', target: document.getElementById('%s'), data: %s, title_color: '%s', header_color: '%s', text_color: '%s', raw: %s"%(width, target, data, title_color, header_color, text_color, self._bool(raw))
        src = """
			(function(){
				Retina.add_renderer({ name: 'paragraph', resource: '""" + self.renderer_resource + """', filename: 'renderer.paragraph.js' });
				Retina.load_renderer('paragraph').then( function () { Retina.Renderer.create('paragraph', {""" + opt + """} ).render(); } );
			})();
		"""
        if self.debug:
            print src
        else:
            IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))
    
    def table(self, width=700, height=600, target="", data=None, rows_per_page=20, sort_autodetect=True, filter_autodetect=True):
        """Table Renderer

          Displays a browsable, filterable table with clickable cells / rows.

          Options

          target (HTML Container Element)
              Element to render the table in.

          width (INT)
              Width of the table.

          height (INT)
              Height of the table.

          rows_per_page (INT)
              The maximum number of table rows to be displayed at a time. Default is 10.

          sortcol (INT)
              Zero based index of the row the table should be sorted by. Default is 0.

          sorted (BOOLEAN)
              Enables / disabled initial sorting of the table by the sortcol. Default is false.

          offset (INT)
              Initial first row to display. Default is 0.

          invisible_columns (HASH)
              Hash of column indices pointing at 1. Columns in this hash are not displayed.

          disable_sort (HASH)
              Hash of column indices pointing at 1. Columns in this hash can not be sorted.

          sorttype (HASH)
              Hash of column indices pointing at a sorttype. A sorttype can be either string or number.

          filter_autodetect (BOOLEAN)
              If set to false will try to detect which filter type is most appropriate for each column. Default is false.

          filter_autodetect_select_max (INT)
              Maximum number of distinct entries in a column that will still autodetec the column filter as a select box. Default is 10.

          sort_autodetect (BOOLEAN)
              If set to true will try to detect which sorttype is appropriate for each column. Default is false.

          filter (HASH)
              Hash of column indices pointing at filter objects. A filter object has the properties
                searchword - the current entry in the search field
                case_sensitive - boolean to turn on / off case sensitivity in filtering
                operator - list of operators available in this filter
                active_operator - selected operator
                type - text or select

          hide_options (BOOLEAN)
              Turns display of the options button on and off. Default is false (the option button is visible).

          onclick (FUNCTION)
              The function to be called when the table is clicked. This function will be passed the parameters (as an ordered list)
                clicked_row - array of contents of the cells of the clicked row
                clicked_cell - content of the clicked cell
                clicked_row_index - zero based index of the clicked row
                clicked_cell_index - zero based index of the clicked cell"""
        if not target:
            target = 'div_'+ipyTools.random_str()
        html = "<div id='%s'></div>"%(target)
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
        if data is None:
            data = "Retina.RendererInstances.table[0].exampleData()"
        else:
            data = json.dumps(data)
            
        opt = "width: %d, height: %d, target: document.getElementById('%s'), data: %s, rows_per_page: %d, sort_autodetect: %s, filter_autodetect: %s"%(width, height, target, data, rows_per_page, self._bool(sort_autodetect), self._bool(filter_autodetect))
        src = """
            (function(){
                Retina.add_renderer({ name: 'table', resource: '""" + self.renderer_resource + """', filename: 'renderer.table.js' });
                Retina.load_renderer('table').then( function () { Retina.Renderer.create('table', {""" + opt + """} ).render(); } );
            })();
        """
        if self.debug:
            print src
        else:
            IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))
    
    def heatmap(self, width=700, height=600, target="", data=None, tree_height=50, tree_width=50, legend_height=250, legend_width=250, row_text_size=15, col_text_size=15, min_cell_height=19):
        """Heatmap Renderer

          Displays a heatmap.

          Options

          width (int)
             Number of pixels the resulting image is wide. Default is 700

          height (int)
             Number of pixels the resulting image is high. This will be adjusted to at least rows * min_cell_height + legend and tree heights. Default is 600.

          tree_height (int)
             Number of pixels the dendogram tree for the columns is high. Default is 50.

          tree_width (int)
             Number of pixels the dendogram tree for the rows is wide. Default is 50.

          legend_height (int)
             Number of pixels for the column names. Default is 250.

          legend_width (int)
             Number of pixels for the row names. Default is 250.

          row_text_size (int)
             Number of pixels of the row text font size. Default is 15.

          col_text_size (int)
             Number of pixels of the column text font size. Default is 15.

          min_cell_height (int)
             Minimum number of pixels a row is high. This may cause the defined height of the resulting image to be overwritten. Default is 19.

          selectedRows (array of boolean)
             Returns an array that has a value of true for all row indices that are currently selected.

          data (object)
             columns (array of string)
                names of the columns
             rows (array of string)
                names of the rows
             colindex (array of int)
                1 based indices of the original column order. This converts the original order (columns) into the one ordered by distance.
             rowindex (array of int)
                1 based indices of the original row order. This converts the original order (rows) into the one ordered by distance.
             coldend (array of array of float)
                distance matrix for the columns
             rowdend
                distance matrix for the rows
             data (array of array of float)
                normalized value matrix"""
        if not target:
            target = 'div_'+ipyTools.random_str()
        html = "<div id='%s'></div>"%(target)
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
        if data is None:
            data = "Retina.RendererInstances.paragraph[0].exampleData()"
        else:
            data = json.dumps(data)
            
        opt = "width: %d, height: %d, target: document.getElementById('%s'), data: %s, tree_height: %d, tree_width: %d, legend_height: %d, legend_width: %d, row_text_size: %d, col_text_size: %d, min_cell_height: %d"%(width, height, target, data, tree_height, tree_width, legend_height, legend_width, row_text_size, col_text_size, min_cell_height)
        src = """
            (function(){
                Retina.add_renderer({ name: 'heatmap', resource: '""" + self.renderer_resource + """', filename: 'renderer.heatmap.js' });
                Retina.load_renderer('heatmap').then( function () { Retina.Renderer.create('heatmap', {""" + opt + """} ).render(); } );
            })();
        """
        if self.debug:
            print src
        else:
            IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))
        
    def _bool(self, aBool):
        if aBool:
            return 'true'
        else:
            return 'false'
