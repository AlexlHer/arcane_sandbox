# Arcane Sandbox
## Personal Arcane Sandbox to test anything

Use main `CMakeLists.txt` to enable or disable SayHello projects.

### Build and test Arcane and Arcane_Sandbox (with pzc: https://github.com/AlexlHer/zsh_config)
```sh
clonearc
initarc D
configarc
biarc

cdw
git clone https://github.com/AlexlHer/arcane_sandbox
initap arcane_sandbox D
ninja

${AP_BUILD_DIR}/template/HelloWorld -A,MaxIteration=1 ${AP_SOURCE_DIR}/template/HelloWorld.arc
${AP_BUILD_DIR}/runcommand/RunCommand -A,MaxIteration=1 ${AP_SOURCE_DIR}/runcommand/HelloWorld.arc
${AP_BUILD_DIR}/itemdirectionmng/ItemDirectionMng -A,MaxIteration=1 ${AP_SOURCE_DIR}/itemdirectionmng/HelloWorld.arc
```
