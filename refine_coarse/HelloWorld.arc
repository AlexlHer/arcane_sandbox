<?xml version="1.0"?>
<case codename="HelloWorld" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3steps</title>
    <timeloop>HelloWorldLoop</timeloop>
  </arcane>


  <mesh amr="true">
    <meshgenerator>
      <cartesian>
        <nsd>2 2</nsd>
        <origine>0.0 0.0</origine>
        <lx nx='32'>32.0</lx>
        <ly ny='32'>32.0</ly>
      </cartesian>
    </meshgenerator>
  </mesh>

  <arcane-checkpoint>
    <checkpoint-service name="ArcaneBasic2CheckpointWriter" />
    <do-dump-at-end>true</do-dump-at-end>
  </arcane-checkpoint>


  <arcane-post-processing>
    <output-period>1</output-period>
    <output>
      <variable>CellCenterCoord</variable>
      <group>AllCells</group>
      <group>AllNodes</group>
    </output>
  </arcane-post-processing>

  <say-hello>
  </say-hello>
</case>
