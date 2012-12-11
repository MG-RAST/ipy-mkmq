#!/usr/bin/env python

import IPython.core.display

class Retina:
	def __init__(self, action='none'):
		self.action = action
		self.renderer_resource = "http://raw.github.com/MG-RAST/Retina/master/renderers/";
		src = """
			(function(){
				Retina.init( { library_resource: "http://raw.github.com/MG-RAST/Retina/master/" });
			})();
		"""
		IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src));
    
	def barchart(self, width='800', height='400', btype="column", target="target", data="Retina.Renderer.graph.exampleData()", title="Browser Usage", x_labels="[ '2005', '2006', '2007', '2008' ]", x_title="Year", y_title="Percentage"):
		html = """
			<div id='""" + target + """'></div>
		"""
		IPython.core.display.display_html(IPython.core.display.HTML(data=html))
		src = """
			(function(){
				Retina.add_renderer({"name": "graph", "resource": "http://raw.github.com/MG-RAST/Retina/master/renderers/", "filename": "renderer.graph.js" });
				Retina.load_renderer("graph").then( function () { Retina.Renderer.graph.render( { width: """ + width + """, height: """ + height + """, type: '""" + btype + """', target: document.getElementById('""" + target  + """'), data: """ + data + """, title: '""" + title + """', x_labels: """ + x_labels + """, x_title: '""" + x_title + """', y_title: '""" + y_title + """' }) });
                        })();
		"""
		IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))

	def paragraph(self, data="Retina.Renderer.paragraph.exampleData()", width="span12", style="", target="target"):
		html = """
			<div id='""" + target + """'></div>
		"""
		IPython.core.display.display_html(IPython.core.display.HTML(data=html))
		src = """
			(function(){
				Retina.add_renderer({ name: 'paragraph', resource: '""" + self.renderer_resource + """', filename: 'renderer.paragraph.js' });
				Retina.load_renderer('paragraph').then( function () { Retina.Renderer.paragraph.render( { target: document.getElementById('""" + target + """'), data: """ + data + """, width: '""" + width + """', style: '""" + style + """' } ) } );
			})();
		"""
		IPython.core.display.display_javascript(IPython.core.display.Javascript(data=src))
