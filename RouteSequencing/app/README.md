## Procedimiento para generar una nueva app

Supongamos que estan instalados Docker, git y tambien el cliente `rc-cli`.

1. Ir a la carpeta donde quiero generar la nueva imagen, por ejemplo, `RC-2021/pruebas/` y hago

```bash
$ rc-cli new-app alpha-centauri rc_python
the 'rc_python' template has been created at '/home/mforets/Projects/RC-2021/prueba/alpha-centauri'
```
Eso crea un nuevo directorio `alpha-centauri` con la estructura de archivos.

2. Substituir los archivos `alpha-centauri/model_apply.sh`, `alpha-centauri/model_build.sh` y `alpha-centauri/Dockerfile` con los correspondientes en la carpeta `RouteSequencing/app`.

3. En la carpeta `alpha-centauri/src` substituir los archivos `model_apply.py` y `model_build.py` por los correspondientes en la carpeta `RouteSequencing/app`. Tambien copiar los archivos `model_apply.jl`, `model_apply.sh`, `model_build.sh`.

4. Copiar en `alpha-centauri/src` el archivo del projecto `Project.toml`, que contiene las dependencias (entre las cuales *no* se incluye `RouteSequencing.jl`) 

5. Copiar en `alpha-centauri/src` la carpeta `RouteSequencing` que contiene el codigo fuente del modulo (archivo `Project.toml` y carpeta `RouteSequencing/src`.

## Probar la app

Una vez completados los pasos anteriores, corremos el Docker (lleva algunos minutos) y calculamos el score.

6. `$ rc-cli model-apply`.

7. `$ rc-cli model-score`.

Para guardar el snapshot: 

8. `$ rc-cli save-snapshot alpha-centauri-v1`.

El production test se hace:

9. `$ rc-cli production-test`
