<?xml version="1.0"?>
<case codename="HelloWorld" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3steps</title>
    <timeloop>HelloWorldLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <save-init>1</save-init>
    <end-execution-output>1</end-execution-output>
    <output-period>1</output-period>
    <output>
      <variable>AMR</variable>
    </output>
  </arcane-post-processing>

  <mesh amr-type="3">
    <meshgenerator>
      <cartesian>
        <nsd>1 1</nsd>
        <origine>0.0 0.0</origine>
        <lx nx='20'>20.0</lx>
        <ly ny='20'>20.0</ly>
      </cartesian>
    </meshgenerator>
  </mesh>

  <say-hello>
  
  </say-hello>
</case>
