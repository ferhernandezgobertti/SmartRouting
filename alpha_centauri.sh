#!/bin/bash
BASE_DIR=$(pwd)
RC_NAME=$1
GIT_CLONE=$2

echo "YOUR CURRENT LOCATION: ${BASE_DIR}"

if [[ $GIT_CLONE -eq "Y" ]]
then
  echo "CLONING REPOSITORY..."
  git clone https://github.com/mforets/RC-2021.git
  echo "REPOSITORY CLONED..."
fi

echo "REPLACING .SH AND DOCKERFILE ..."
cp -f "${BASE_DIR}/RC-2021/RouteSequencing/app/model_apply.sh" "${BASE_DIR}/${RC_NAME}/model_apply.sh" 
cp -f "${BASE_DIR}/RC-2021/RouteSequencing/app/model_build.sh" "${BASE_DIR}/${RC_NAME}/model_build.sh" 
cp -f "${BASE_DIR}/RC-2021/RouteSequencing/app/Dockerfile" "${BASE_DIR}/${RC_NAME}/Dockerfile" 

echo "REPLACING APPLY.PY AND BUILD.PY ..."
cp -f "${BASE_DIR}/RC-2021/RouteSequencing/app/src/model_apply.py" "${BASE_DIR}/${RC_NAME}/src/model_apply.py" 
cp -f "${BASE_DIR}/RC-2021/RouteSequencing/app/src/model_build.py" "${BASE_DIR}/${RC_NAME}/src/model_build.py" 

echo "ACQUIRING APPLY.JL, APPLY.SH AND BUILD.SH ..."
cp "${BASE_DIR}/RC-2021/RouteSequencing/app/src/model_apply.jl" "${BASE_DIR}/${RC_NAME}/src/" 
cp "${BASE_DIR}/RC-2021/RouteSequencing/app/src/model_apply.sh" "${BASE_DIR}/${RC_NAME}/src/" 
cp "${BASE_DIR}/RC-2021/RouteSequencing/app/src/model_build.sh" "${BASE_DIR}/${RC_NAME}/src/" 

echo "ACQUIRING PROJECT.TOML ..."
cp "${BASE_DIR}/RC-2021/RouteSequencing/app/src/Project.toml" "${BASE_DIR}/${RC_NAME}/src/" 

echo "ACQUIRING ROUTESEQUENCING FOLDER ..."
cp -r "${BASE_DIR}/RC-2021/RouteSequencing" "${BASE_DIR}/${RC_NAME}/src/" 

echo "DONE"