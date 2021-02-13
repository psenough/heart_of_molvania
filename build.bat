rmdir /s /q dist

mkdir dist
mkdir dist\gfx
mkdir dist\audio

xcopy /s gfx dist\gfx
xcopy /s audio dist\audio

java -XX:ReservedCodeCacheSize=64m -jar bin/html-compressor.jar --compress-css --compress-js -o dist/index.html index.html
java -XX:ReservedCodeCacheSize=64m -jar bin/js-compiler.jar --js utility.js --js_output_file dist/utility.js

