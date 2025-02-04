<?xml version="1.0"?>
<case codename="HelloWorld" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3steps</title>
    <timeloop>HelloWorldLoop</timeloop>
  </arcane>

  <mesh amr="true">
    <meshgenerator>
      <cartesian>
        <nsd>8 8</nsd>
        <origine>0.0 0.0</origine>
        <lx nx='8'>8.0</lx>
        <ly ny='8'>8.0</ly>
      </cartesian>
    </meshgenerator>
  </mesh>

  <arcane-checkpoint>
  </arcane-checkpoint>
  <say-hello>
  </say-hello>
</case>

<!--
./test_replace/TestReplace \
-A,arcane-checkpoint/checkpoint-service/@name=ArcaneBasic2CheckpointWriter \
-A,arcane-checkpoint/do-dump-at-end=true \
-A,say-hello/test-option=1 \
-A,say-hello/boundary-condition=Y \
-A,say-hello/pdes-random-number-generator/@name=PDESRandomNumberGenerator \
-A,say-hello/pdes-random-number-generator/@mesh-name=Mesh0 \
-A,say-hello/pdes-random-number-generator/initialSeed=123456 \
${AP_SOURCE_DIR}/test_replace/HelloWorld2.arc
-->