#!/usr/bin/env python

import IPython.core.display

class Retina:
    def __init__(self, action='none', debug=False):
        self.action = action
        self.debug  = debug
        self.renderer_resource = "http://raw.github.com/MG-RAST/Retina/master/renderers/";
        src = """
			(function(){
				Retina.init( { library_resource: "http://raw.github.com/MG-RAST/Retina/master/js/" });
			})();
		"""
        IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))

    def graph(self, width=800, height=400, btype="column", target="target", data="Retina.RendererInstances.graph[0].exampleData()", title="Browser Usage", x_labels="['2005','2006','2007','2008']", x_title="Year", y_title="Percentage", show_legend=False, legend_position='left'):
        html = """
			<div id='""" + target + """'></div>
		"""
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
        opt = "width: %d, height: %d, type: '%s', target: document.getElementById('%s'), data: %s, title: '%s', x_labels: %s, x_title: '%s', y_title: '%s', show_legend: %s, legend_position: '%s'"%(width, height, btype, target, data, title, x_labels, x_title, y_title, self._bool(show_legend), legend_position)
        src = """
			(function(){
				Retina.add_renderer({"name": "graph", "resource": '""" + self.renderer_resource + """', "filename": "renderer.graph.js" });
				Retina.load_renderer("graph").then( function () { Retina.Renderer.create('graph', {""" + opt + """}).render(); });
                        })();
		"""
        if self.debug:
            print src
        else:
            IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))
    
    def plot(self, width=800, height=400, target="target", data="Retina.RendererInstances.plot[0].exampleData()", title="Sine", show_legend=True, legend_position='left'):
        html = """
			<div id='""" + target + """'></div>
		"""
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
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

    def paragraph(self, width="span12", target="target", data="Retina.RendererInstances.paragraph[0].exampleData()", title_color='black', header_color='black', text_color='black', raw=False):
        html = """
			<div id='""" + target + """'></div>
		"""
        IPython.core.display.display_html(IPython.core.display.HTML(data=html))
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

    def _bool(self, aBool):
        if aBool:
            return 'true'
        else:
            return 'false'
