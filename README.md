# PolarCodes
Keeps code related to the polar code project




COMPILE MAC: 
/usr/local/bin/cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_C_COMPILER:FILEPATH=/usr/local/bin/gcc-10 -DCMAKE_CXX_COMPILER:FILEPATH=/usr/local/bin/g++-10 -H/Users/valyagorbunova/Desktop/Работа/huawei/PolarCodes -B/Users/valyagorbunova/Desktop/Работа/huawei/PolarCodes/build -G Ninja

remove -fno-fat-lto-objects from build/ninja.build

/usr/local/bin/cmake --build /Users/valyagorbunova/Desktop/Работа/huawei/PolarCodes/build --config Release --target all -- -j 6   

RUN MAC: ./x64/Release/1PolarSimulate /Users/valyagorbunova/Desktop/Работа/huawei/PolarCodes/configs/test_config.xml  