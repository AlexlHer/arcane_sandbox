<?xml version="1.0"?>
<case codename="HelloWorld" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3steps</title>
    <timeloop>HelloWorldLoop</timeloop>
  </arcane>

  <meshes >
    <mesh nb-ghostlayer="0">
      <generator name="Cartesian3D" >
        <nb-part-x>1</nb-part-x> 
        <nb-part-y>2</nb-part-y>
        <nb-part-z>2</nb-part-z>
        <origin>0.0 0.0 0.0</origin>

        <x>
          <n>10</n>
          <length>10</length>
        </x>

        <y>
          <n>2</n>
          <length>2</length>
        </y>

        <z>
          <n>2</n>
          <length>2</length>
        </z>

      </generator>
    </mesh>
  </meshes>

  <say-hello>

  </say-hello>
</case>
