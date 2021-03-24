from aiohttp import web
import socketio
import os
import time


sio = socketio.AsyncServer(cors_allowed_origins="*")
app = web.Application()
sio.attach(app)


async def index(request):
    """Serve the client-side application."""
    return web.Response(text="It's Alive!", content_type="text/html")


@sio.on("get pypka status", namespace="/")
async def pypka_status(sid, subID):

    logpath = f"/tmp/tmp_{subID}/LOG_{subID}"

    content = []
    if os.path.isfile(logpath):
        with open(logpath) as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith("PB Runs"):
                    if content[-1].startswith("PB Runs"):
                        content[-1] = line
                    else:
                        content.append(line)
                elif line.startswith("MC Run"):
                    if content[-1].startswith("MC Run"):
                        content[-1] = line
                    else:
                        content.append(line)
                else:
                    content.append(line)
            content = "\n".join(content)

    await sio.emit(
        "pypka status", {"content": content if content else None, "subID": subID}
    )
    print("finished", subID, len(content))


app.router.add_get("/", index)

if __name__ == "__main__":
    web.run_app(app)
