# CelineEngine

The circuit simulator shared by Céline and Céline Player.

Modified nodal analysis with SPICE-grade device models — diodes, BJTs, JFETs,
MOSFETs, op-amps, triodes, pentodes, rectifiers, transformers — solved per sample
by the DK method: the linear part of the circuit is eliminated once, and Newton runs
only on the ports of the nonlinear devices. A netlist goes in, audio comes out; it
knows nothing about drawing or plugins.

Read the class comment in `CelineEngine/Engine.h` first, and `LIMITATIONS.md` before
trusting a result.

## Using it from a plugin

As a submodule, beside JUCE and the others:

```sh
git submodule add https://github.com/Celine-audio/CelineEngine.git modules/CelineEngine
```

In the plugin's `CMakeLists.txt`, once `SharedCode` exists:

```cmake
set(CELINE_ENGINE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/modules/CelineEngine" CACHE PATH "The CelineEngine checkout to build against")
include("${CELINE_ENGINE_DIR}/CelineEngine.cmake")
target_link_libraries(SharedCode INTERFACE CelineEngine)
```

Headers are included by repository: `#include <CelineEngine/Engine.h>`. It needs
JUCE's `juce_audio_processors` and nothing from the plugin.

```cpp
Circuit circuit;
circuit.addResistor ("in", "out", 1000.0);
circuit.addCapacitor ("out", "gnd", 100e-9);
circuit.setInputNode ("in");
circuit.setOutputNode ("out");
circuit.prepare (48000.0);

float y = circuit.process (x);
```

## Working on it

Edit a checkout beside the plugins and build one against it with
`-DCELINE_ENGINE_DIR=../CelineEngine`. Commit and push here first, then move each
plugin's submodule to that commit. The tests are in Céline, which can also build
the drawn sheets that exercise the engine end to end.

## Licence

AGPLv3, like the plugins that compile it. See `LICENSE` and `COPYING`.
