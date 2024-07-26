<?xml version="1.0"?>
<case codename="HelloWorld" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3steps</title>
    <timeloop>HelloWorldLoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <generator name="Cartesian2D" >
        <nb-part-x>1</nb-part-x> 
        <nb-part-y>1</nb-part-y>
        <origin>0.0 0.0</origin>
        <x><n>20</n><length>2.0</length></x>
        <y><n>20</n><length>2.0</length></y>
      </generator>
    </mesh>
  </meshes>

  <say-hello>
    <nb-sds>100</nb-sds>
    <nb-blocs>20</nb-blocs>
    <nb-cells>30</nb-cells>
    <nb-nodes>100</nb-nodes>
    <nb-values>100</nb-values>
  </say-hello>
</case>
