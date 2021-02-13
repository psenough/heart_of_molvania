JAVA_FLAGS    = -XX:ReservedCodeCacheSize=64m
JS_COMPILER   = 'bin/js-compiler.jar'
HTML_COMPILER = 'bin/html-compressor.jar'

HTMLS = 							\
	./index.html
HTMLS_MIN = $(subst ./,dist/,${HTMLS})

SCRIPTS = 						    \
	./utility.js
SCRIPTS_MIN = $(subst ./,dist/,${SCRIPTS})

SHADERS = 						    \
	./frag.glsl
SHADERS_MIN = $(subst ./,dist/,${SHADERS})

all: html js res shaders

clean:
	@@rm -rf dist

init:
	@@mkdir -p dist

res:
	@@mkdir -p dist/audio
	@@mkdir -p dist/gfx
	@@cp gfx/* dist/gfx
	@@cp audio/* dist/audio

dist/%.html : %.html
	@@java ${JAVA_FLAGS} -jar ${HTML_COMPILER} --compress-css --compress-js -o $@ $<

dist/%.js : %.js
	@@java ${JAVA_FLAGS} -jar ${JS_COMPILER} --js $< --js_output_file $@

dist/%.glsl : %.glsl
	@@cp $< $@


html: init $(HTMLS_MIN)

js: init $(SCRIPTS_MIN)

shaders: init $(SHADERS_MIN)


