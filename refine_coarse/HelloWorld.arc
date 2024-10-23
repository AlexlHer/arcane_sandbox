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
    <checkpoint-service name="ArcaneBasic2CheckpointWriter" />
    <do-dump-at-end>true</do-dump-at-end>
  </arcane-checkpoint>

  <say-hello>
  
  </say-hello>
</case>
