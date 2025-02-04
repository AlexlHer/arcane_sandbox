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
    <checkpoint-service name="@test5@" />
    <do-dump-at-end>true</do-dump-at-end>
  </arcane-checkpoint>

  <say-hello>
    <test-option>@test@</test-option>
    <boundary-condition>@test2@</boundary-condition>
    <pdes-random-number-generator name="@test3@" mesh-name="@mesh_name@">
      <initialSeed>@test4@</initialSeed>
    </pdes-random-number-generator>
  </say-hello>
</case>

<!--
./test_replace/TestReplace \
-A,test5=ArcaneBasic2CheckpointWriter \
-A,test2=X \
-A,test3=PDESRandomNumberGenerator \
-A,test4=678 \
-A,test6=ArcaneBasic2CheckpointWriter \
-A,say-hello/test-option=@test4@ \
-A,say-hello/pdes-random-number-generator/@mesh-name=Mesh0 \
${AP_SOURCE_DIR}/test_replace/HelloWorld.arc
-->