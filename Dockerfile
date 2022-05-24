FROM eclipse-temurin:11-alpine

RUN apk update && \
        apk upgrade && \
        apk add --no-cache \
            cmake=3.21.3-r0 \
            make \
            g++ \
            gcc \
            git

ADD /HiGHS.zip /
RUN unzip HiGHS.zip -d /HiGHS
RUN mkdir HiGHS/build && \
        cd HiGHS/build && \
        cmake .. && \
        cmake --build .

# folders should have the 775 permission, which is the minimum allowed permissions to write the files within the API
RUN mkdir HiGHS/build/bin/input && \
        mkdir HiGHS/build/bin/output && \
        chmod 775 HiGHS/build/bin/input && \
        chmod 775 HiGHS/build/bin/output
