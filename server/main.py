import asyncio
from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from server.middlewares.exception_handlers import catch_exception_middeware

async def initialize_deseq2_service():
    from server.routes.run_deseq2 import router as run_deseq2

    app.include_router(run_deseq2)


@asynccontextmanager
async def lifespan(app: FastAPI):
    asyncio.create_task(initialize_deseq2_service())
    yield

app = FastAPI(title='Deseq2 Runner', description='Run DeSeq2 Analysis', lifespan=lifespan)

# CORS Setup
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], 
    allow_credentials=["*"], 
    allow_methods=["*"], 
    allow_headers=["*"]
)

# Middleware exception handlers
app.middleware("http")(catch_exception_middeware)